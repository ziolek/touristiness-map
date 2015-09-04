#!/usr/bin/python

# This program creates the "Touristiness" and "Interesting remote places"
# map overlay images. For its input, it requires
#	* GeoNames cities datafile, obtained from
#			http://download.geonames.org/export/dump/cities1000.zip
#			(unzip the file before use)
#	* Internet access for queries to Panoramio API
#
# Panoramio API query results are optionally cached on disk between runs;
# this requires ca 340Mbytes of disk space. If cache is empty or disabled,
# the program runs for ca 4 hours (Panoramio requests are done in parallel
# with Python threads); with cache filled (if all Panoramio queries are
# already done) it runs for ca 2 minutes. RAM usage is ca 100Mbytes.
#
# Output images are 720x720 paletted PNGs with alpha channel in palette.
# Resolution is 1/4 coordinate degrees, which results in 1440x720 pixels
# for the whole world. Since Google Maps does not support overlay images
# spanning the whole globe, the world is broken down into two 720x720
# images, one for each hemisphere.

from __future__ import with_statement
import sys,math,operator,re,time,cPickle,threading,socket,httplib,libxml2
import png		# from http://packages.python.org/pypng/

mintreelevel=2
maxtreelevel=5

if len(sys.argv) < 1+4:
	print 'Usage: build_map.py <cachedir (empty string disables caching)> ' + \
			'<geonames cities1000.txt input filename> ' + \
			'<touristiness map output filename prefix> ' + \
			'<interesting remote places map output filename prefix>'
	sys.exit(1)

cachedir=sys.argv[1]
geonames_cities_fname=sys.argv[2]
touristiness_image_fname_prefix=sys.argv[3]
remote_interesting_image_fname_prefix=sys.argv[4]
print_debug=False

thread_local_data=threading.local()
threads_stop_flag=False

touristy_spots=[]
touristy_spots_lock=threading.Lock()

slots_per_degree=(2**mintreelevel)

########################################################################
######################## Panoramio API access ##########################
########################################################################

libxml2.registerErrorHandler(lambda x,y:x,None)

panoramio_hostname='www.panoramio.com'
panoramio_ipaddr=socket.gethostbyname(panoramio_hostname)

def do_http_request(url):
	global thread_local_data,threads_stop_flag
	global panoramio_hostname,panoramio_ipaddr

	for try_nr in range(10):
		if threads_stop_flag:
			return None
		try:
			if not hasattr(thread_local_data,'conn'):
				thread_local_data.conn=httplib.HTTPConnection(
														panoramio_ipaddr)
			thread_local_data.conn.request('GET',url,
									headers={'Host':panoramio_hostname})
			r=thread_local_data.conn.getresponse()
			data=r.read()
			if r.status == 200:
				break
			if (r.status >= 400 and r.status < 500) or r.status in (301,302):
				print >> sys.stderr,url,'received HTTP error',r.status,r.reason
				return None
			print >> sys.stderr,url,'Received HTTP error',r.status,r.reason
		except (httplib.BadStatusLine,httplib.CannotSendRequest,
										ValueError,socket.gaierror),err:
			print >> sys.stderr,url,'httplib error:',err
		thread_local_data.conn.close()
		del thread_local_data.conn
		print >> sys.stderr,url,'Pausing...'
		time.sleep(20)
	else:
		return None

	data=re.sub(r'<django.templatetags.i18n.TranslateNode [^>]*>','',data)
	try:
		return libxml2.parseDoc(data)
	except libxml2.parserError:
		print >> sys.stderr,'Error parsing this XML response:\n',data
		return None

username_re=re.compile('panoramio.com/user/([0-9]+)')

def fetch_nr_of_photos(lat_range,lon_range):
	global print_debug
	doc=do_http_request('/panoramio.kml?BBOX=%s,%s,%s,%s' % \
					(lon_range[0],lat_range[0],lon_range[1],lat_range[1]))
	if print_debug:
		sys.stdout.write('.')
		sys.stdout.flush()

	if not doc:
		return None
	ctx=doc.xpathNewContext()
	ctx.xpathRegisterNs('n',doc.children.ns().content)
	userids=[]
	for placemark in ctx.xpathEval('/n:kml/n:Document/n:Placemark'):
		ctx1=doc.xpathNewContext()
		ctx1.setContextNode(placemark)
		ctx1.xpathRegisterNs('n',doc.children.ns().content)
		description=' '.join(map(str,ctx1.xpathEval('n:description/text()')))
		ctx1.xpathFreeContext()
		r=username_re.search(description)
		if r:
			userids.append(int(r.group(1)))
		else:
			userids.append(hash(str(placemark)))
	ctx.xpathFreeContext()
	doc.freeDoc()
	return userids

def coords_to_idxs(lat,lon):
	global slots_per_degree
	return (int(math.floor((lat - (-90 )) * slots_per_degree)),
			int(math.floor((lon - (-180)) * slots_per_degree)))

########################################################################
########################## outimage_writer #############################
########################################################################

outimages=[]

class outimage_writer:
	def __init__(self,lat_range,lon_range,fname_prefix,palette):
		global slots_per_degree,outimages
		self.lat_range,self.lon_range=lat_range,lon_range
		self.start_idxs=coords_to_idxs(self.lat_range[0],self.lon_range[0])
		self.size=(
				(self.lon_range[1]-self.lon_range[0]) * slots_per_degree,
				(self.lat_range[1]-self.lat_range[0]) * slots_per_degree)
		self.imagedata=[0] * (self.size[0]*self.size[1])
		self.palette=palette
		self.fname='%s_%d_%d_%d_%d.png' % (fname_prefix,
									self.lat_range[0],
									self.lon_range[0],
									self.lat_range[1]-self.lat_range[0],
									self.lon_range[1]-self.lon_range[0])
		self.lock=threading.Lock()
		outimages.append(self)

	def latlon_to_outimage_coords(self,lat,lon):
		global slots_per_degree
		return (int(round((lon-self.lon_range[0]) * slots_per_degree)),
				int(round((lat-self.lat_range[0]) * slots_per_degree)))

	def putpixel(self,x,y,coloridx):
		if x >= 0 and y >= 0 and x < self.size[0] and y <= self.size[1]:
			with self.lock:
				self.imagedata[x + (self.size[1]-1-y)*self.size[0]]=coloridx

	def putrectangle(self,x0,x1,y0,y1,coloridx):
		with self.lock:
			for x in range(max(0,x0),min(self.size[0],x1)):
				for y in range(max(0,self.size[1]-1-y1),min(self.size[1],self.size[1]-1-y0)):
					self.imagedata[x + y*self.size[0]]=coloridx

	def write_file(self):
		f=open(self.fname,'wb')
		png.Writer(size=self.size,palette=self.palette,compression=9). \
											write_array(f,self.imagedata)
		f.close()

########################################################################
######################### Population map ###############################
########################################################################

population_map_horizontal_sum=[]
population_map_vertical_sum=[]

def build_population_map_sums():
	global population_map,slots_per_degree
	global population_map_vertical_sum,population_map_horizontal_sum

	population_map_vertical_sum[:]=[]
	cur_sums=[0] * 360 * slots_per_degree
	for lat in range(len(population_map)):
		for lon in range(len(cur_sums)):
			cur_sums[lon]+=population_map[lat][lon]
		population_map_vertical_sum.append(tuple(cur_sums))
	population_map_vertical_sum.append(tuple([0] * 360 * slots_per_degree))

	population_map_horizontal_sum[:]=[[0]*len(population_map[0]) \
								for lat in range(len(population_map))]
	cur_sums=[0] * 180 * slots_per_degree
	for lon in range(len(population_map[0])):
		for lat in range(len(cur_sums)):
			cur_sums[lat]+=population_map[lat][lon]
			population_map_horizontal_sum[lat][lon]=cur_sums[lat]
	for lat in range(len(cur_sums)):
		population_map_horizontal_sum[lat].append(0)

def build_empty_map():
	global slots_per_degree
	m=[]
	for lat in range(-90 * slots_per_degree,+90 * slots_per_degree):
		m.append([0] * 360 * slots_per_degree)
	return m

print 'Reading GeoNames file'

population_map=build_empty_map()

for line in open(geonames_cities_fname,'r'):
	fields=line.split('\t')
	lat,lon=map(float,fields[4:6])
	population=float(fields[14])

	lat,lon=coords_to_idxs(lat,lon)
	population_map[lat][lon]+=int(population)

build_population_map_sums()

def calc_remoteness(lat,lon):
	global slots_per_degree,population_map
	global population_map_horizontal_sum,population_map_vertical_sum

	max_distance_from_city=4 * slots_per_degree
	best_distance_from_city=max_distance_from_city
	distance=0
	y_size=len(population_map[0])
	while distance < best_distance_from_city+1*slots_per_degree:
		if distance:

			if False:
				population=0
				x=lat-distance
				if x >= 0:
					for y in range(lon-distance,lon+distance+1):
						population+=population_map[x][y % y_size]
				x=lat+distance
				if x < 180 * slots_per_degree:
					for y in range(lon-distance,lon+distance+1):
						population+=population_map[x][y % y_size]
				y1=(lon-distance) % y_size
				y2=(lon+distance) % y_size
				for x in range(max(0,lat-distance+1),
								min(180 * slots_per_degree,lat+distance)):
					population+=population_map[x][y1]+population_map[x][y2]
			else:
				y1=(lon-distance) % y_size
				y2=(lon+distance) % y_size
				x1=max(-1,lat-distance)
				x2=min(180 * slots_per_degree-1,lat+distance-1)
				population=	(population_map_vertical_sum[x2][y2] - \
							population_map_vertical_sum[x1][y2]) + \
							(population_map_vertical_sum[x2][y1] - \
							population_map_vertical_sum[x1][y1])

				y1-=1
				if y2 > y1:
					y3,y4=0,0
				else:
					y3,y4=-1,y2
					y2=y_size-1
	
				x=lat-distance
				if x >= 0:
					population+=(population_map_horizontal_sum[x][y2] - \
								population_map_horizontal_sum[x][y1]) + \
								(population_map_horizontal_sum[x][y4] - \
								population_map_horizontal_sum[x][y3])

				x=lat+distance
				if x < 180 * slots_per_degree:
					population+=(population_map_horizontal_sum[x][y2] - \
								population_map_horizontal_sum[x][y1]) + \
								(population_map_horizontal_sum[x][y4] - \
								population_map_horizontal_sum[x][y3])
		else:
			population=population_map[lat][lon]

		if population:
			distance_from_city=distance
			if population >= 50e3:
				distance_from_city-=slots_per_degree * \
								math.sqrt(max(0,population - 50e3) / 60e6)
			else:
				distance_from_city+=slots_per_degree * 2000.0 / population
			best_distance_from_city=min(best_distance_from_city,
														distance_from_city)
		distance+=1

	if best_distance_from_city <= 0:
		return 0
	return min(1,best_distance_from_city / float(max_distance_from_city))

########################################################################
####################### toplevel_rectangle_processor ###################
########################################################################

class toplevel_rectangle_processor:
	def __init__(self,lat,lon,touristiness_outimage):
		global cachedir
		self.lat=lat
		self.lon=lon
		self.touristiness_outimage=touristiness_outimage

		self.cachefname=(cachedir + '/%d_%d.pickle' % (lat,lon)) \
													if cachedir else None
		self.cache=dict()
		if self.cachefname:
			try:
				loaded_cache=cPickle.load(open(self.cachefname,'r'))
			except (IOError,EOFError):
				pass
			else:
				self.cache=loaded_cache
		self.cache_modifications=0

	def get_photo_userids(self,lat_range,lon_range):
		global print_debug

		photo_userids=self.cache.get((lat_range,lon_range))
		if photo_userids is None:
			photo_userids=fetch_nr_of_photos(lat_range,lon_range)
			if photo_userids is not None:
				self.cache[lat_range,lon_range]=photo_userids
				self.cache_modifications+=1
				self.check_write_cache(50)
		elif print_debug:
			sys.stdout.write('c')
			sys.stdout.flush()
		return photo_userids or tuple()

	def check_write_cache(self,min_modifications=1):
		if self.cache_modifications >= min_modifications and self.cachefname:
			cPickle.dump(self.cache,open(self.cachefname,'w'),-1)
			self.cache_modifications=0

	@staticmethod
	def write_output_rectangle(outimage,lat_range,lon_range,coloridx):
		x0,y0=outimage.latlon_to_outimage_coords(lat_range[0],lon_range[0])
		x1,y1=outimage.latlon_to_outimage_coords(lat_range[1],lon_range[1])
		outimage.putrectangle(x0,x1,y0,y1,coloridx)

	def process_rectangle(self,lat,lon,treelevel=0):
		global touristy_spots,touristy_spots_lock

		coord_step=0.5 ** treelevel
		lat_range=(lat,lat+coord_step)
		lon_range=(lon,lon+coord_step)

		photo_userids=self.get_photo_userids(lat_range,lon_range)

		if photo_userids and treelevel < maxtreelevel and \
					(len(photo_userids) >= 25 or treelevel < mintreelevel):
			nr_of_photos=[]
			for lat_idx in (0,1):
				for lon_idx in (0,1):
					nr_of_photos.append(self.process_rectangle(
									lat+lat_idx*coord_step*0.5,
									lon+lon_idx*coord_step*0.5,treelevel+1))
			nr_of_photos=sum(sorted(nr_of_photos,reverse=True)[:2]) * 2
		else:
			nr_of_userids=len(frozenset(photo_userids))
			nr_of_photos=nr_of_userids + \
								(len(photo_userids) - nr_of_userids) / 3.0

		if not photo_userids and treelevel < mintreelevel:
			self.write_output_rectangle(self.touristiness_outimage,
												lat_range,lon_range,0)
		elif treelevel == mintreelevel:
			level=0
			if nr_of_photos > 0:
				area=math.cos(math.radians(lat+coord_step/2.0)) * \
						((2**2) ** maxtreelevel) / ((2**2) ** treelevel)
				touristiness=nr_of_photos / (30.0 * area)
				with touristy_spots_lock:
					touristy_spots.append([coords_to_idxs(lat,lon),
															touristiness])
				level=int(min(99,max(1,
								99*(1+math.log(min(1,touristiness))/7))))

			x,y=self.touristiness_outimage.latlon_to_outimage_coords(
												lat_range[0],lon_range[0])
			self.touristiness_outimage.putpixel(x,y,level)

		return nr_of_photos

	def process(self):
		self.process_rectangle(self.lat,self.lon)
		self.check_write_cache()

class toplevel_rectangle_processor_thread(threading.Thread,
											toplevel_rectangle_processor):
	def __init__(self,*args):
		threading.Thread.__init__(self)
		toplevel_rectangle_processor.__init__(self,*args)

	def run(self):
		self.process()

########################################################################
###################### Build dest regions list #########################
########################################################################

regions=[
		#((-70,70),(-180,0)),
		#((-60,+33),(-118,-33)),
		#((-6,53),(72,115)),
		]

lat_grid=(-90,+90,180)
lon_grid=(-180,+180,180)

for lat_start in range(*lat_grid):
	for lon_start in range(*lon_grid):
		regions.append(((lat_start,lat_start+lat_grid[2]),
						(lon_start,lon_start+lon_grid[2])))

########################################################################
######################## Build touristiness map ########################
########################################################################

print 'Building touristiness map'

touristiness_color_gradient=[(0,0,0,160)]
for level100 in range(1,100):
	level=level100/(100.0-1)
	touristiness_color_gradient.append(tuple(map(int,(
											min(255,level*2*255),	# R
											max(0,  level*2-1)*255,	# G
											max(0,  1-level*2)*255,	# B
											160+(228-160)*level))))	# A

try:
	for lat_range,lon_range in regions:
		touristiness_outimage=outimage_writer(lat_range,lon_range,
									touristiness_image_fname_prefix,
									touristiness_color_gradient)
		for lon in range(*lon_range):
			for lat in range(*lat_range):
				while threading.activeCount() > 50:
					time.sleep(0.3)
				if print_debug:
					print lon,lat,
					sys.stdout.flush()
				processor_args=(lat,lon,touristiness_outimage)
				if True:
					thread=toplevel_rectangle_processor_thread(
															*processor_args)
					thread.start()
				else:
					toplevel_rectangle_processor(*processor_args).process()
				if print_debug:
					print
	print 'Waiting for fetching threads to finish'
	while threading.activeCount() > 1:
		time.sleep(0.3)
	while outimages:
		outimages.pop().write_file()
except:
	threads_stop_flag=True
	raise

########################################################################
################## Build remote interesting places map #################
########################################################################

print 'Building remote interesting places map'

for idx in xrange(len(touristy_spots)):
	(lat_idx,lon_idx),touristiness=touristy_spots[idx]
	remoteness=calc_remoteness(lat_idx,lon_idx)
	touristy_spots[idx].append(remoteness)
	touristy_spots[idx].append(200e3 * (max(0,touristiness-0.20)**2) \
												if remoteness > 0.1 else 0)

for idx in xrange(len(touristy_spots)):
	(lat_idx,lon_idx),touristiness,remoteness,tourist_population= \
													touristy_spots[idx]
	population_map[lat_idx][lon_idx]+=tourist_population

build_population_map_sums()

interestingness_map=build_empty_map()
nontouristy_interestingness_map=build_empty_map()

filter_coeffs=[]
for lat_idx in range(len(interestingness_map)):
	coeffs=[]
	filter_radius=2
	filter_radius_y=1 + ((filter_radius-1) / \
			math.cos(math.radians((lat_idx + 0.5)/slots_per_degree - 90)))
	for x in range(lat_idx-filter_radius+1,lat_idx+filter_radius):
		if x < 0 or x >= len(interestingness_map):
			continue
		x_partvalue=((x-lat_idx)/float(filter_radius)) ** 2
		for y_delta in range(-int(filter_radius_y),int(filter_radius_y)+1):
			coeff=1 - math.sqrt(x_partvalue + \
								(y_delta / float(filter_radius_y)) ** 2)
			if coeff >= 0.03:
				coeffs.append([x,y_delta,coeff])

	coeffs_sum=sum(map(operator.itemgetter(2),coeffs))
	for idx in range(len(coeffs)):
		coeffs[idx][2]/=coeffs_sum

	filter_coeffs.append(coeffs)

remote_interesting_palette=[]
for i in range(16):
	for j in range(16):
		remote_interesting_palette.append((	int(round(255*i/15)),
											int(round(255*j/15)),
											0,100+80*max(i,j)/16))

for lat_range,lon_range in regions:
	outimage_writer(lat_range,lon_range,
								remote_interesting_image_fname_prefix,
								remote_interesting_palette)

for (lat_idx,lon_idx),touristiness,remoteness,tourist_population in \
															touristy_spots:
	if remoteness <= 0.075:
		continue
	remoteness_withtourists=calc_remoteness(lat_idx,lon_idx)

	#level=int(max(0,min(99,remoteness_withtourists*200)))
	#level=int(min(99,max(1,99*(1+math.log(min(1,touristiness))/7))))
	#for image in outimages:
	#	image.putpixel(	lon_idx-image.start_idxs[1],
	#					lat_idx-image.start_idxs[0],
	#					level)
	#continue

	nontouristy_interestingness=min(1,remoteness_withtourists*3) * \
						(remoteness_withtourists/remoteness + \
						2*min(1,-math.log(min(1,touristiness)) / 7)) / 3.0

	norm_touristiness=1-0.8*min(1,nontouristy_interestingness*2)
	interestingness=min(1,touristiness/0.005) * (1 + math.log(min(1,
				min(1,touristiness/norm_touristiness) * \
				(remoteness**3) * 1000)) / 7)

	nontouristy_interestingness*=interestingness

	for x,y_delta,filter_coeff in filter_coeffs[lat_idx]:
		y=(lon_idx+y_delta) % len(interestingness_map[0])
		interestingness_map[x][y]+=interestingness*filter_coeff
		nontouristy_interestingness_map[x][y]+= \
								nontouristy_interestingness*filter_coeff

for lat_idx in range(len(interestingness_map)):
	for lon_idx,(interestingness,nontouristy_interestingness) in \
							enumerate(zip(  interestingness_map[lat_idx],
								nontouristy_interestingness_map[lat_idx])):
		if interestingness <= 1e-3:
			continue

		brightness=15*min(1,2*interestingness)
		angle=min(1,max(0,nontouristy_interestingness / interestingness))
		coloridx=16*int(round(min(1,2-2*angle) * brightness)) + \
					int(round(min(1,  2*angle) * brightness))

		for image in outimages:
			image.putpixel(	lon_idx-image.start_idxs[1],
							lat_idx-image.start_idxs[0],
							coloridx)

while outimages:
	outimages.pop().write_file()
