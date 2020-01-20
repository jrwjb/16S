#!/share/apps/library/bin/python
from __future__ import print_function,division
from collections import defaultdict
from optparse import OptionParser
from math import sin,cos,pi
from getline_1_1 import levelline
import re,sys,os

p = OptionParser()
p.add_option('-f',dest='file',help="normally , it is otu_taxa_table.xls")
p.add_option('-t',dest='top',help="specify the top number of species to draw if not set '-p'",type="int")
p.add_option('-p',dest='percent',help="specify percent of species above to draw if not set '-t'",type="float")
p.add_option('-o',dest='out',help='the path of out')
p.add_option('--config',dest='config',help="specify the height of each sample, default set 3000 for each sample> 10")
(options,args) = p.parse_args()


def tree(): return defaultdict(tree)

file = open(options.file)
taxline = tree()
phylum = defaultdict(int)
clas =  defaultdict(int)
order = defaultdict(int)
family = defaultdict(int)
genus = defaultdict(int)


samindex = defaultdict(int)

samindex = {}
samkindom = {}
samphylum =	{}
samclas = {}
samorder = {}
samfamily = {}
samgenus = {}
# samdomain = {}
samtaxstop = {}

seqs = 0
title = re.split(r'\t',file.readline().strip())[1:-1]

for i in title:
	samindex[i] = title.index(i)

	samkindom[i] = defaultdict(int)
	samphylum[i] = defaultdict(int)
	samclas[i] = defaultdict(int)
	samorder[i] = defaultdict(int)
	samfamily[i] = defaultdict(int)
	samgenus[i] = defaultdict(int)
	# samdomain[i] = defaultdict(int)


	samtaxstop[i] = defaultdict(int)

samall=defaultdict(int)
samphylumper={}
samorderper={}
samfamilyper={}
samclasper={}
samgenusper={}
domain=''
for i in file.readlines():
	tmp = levelline(i.strip())
	if tmp:
		taxline[tmp[0]][tmp[1]][tmp[2]][tmp[3]][tmp[4]] = 1
		linedata = [int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]]
		# m = re.match(r'.*(d__(?P<domain>\w+))(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())
		m = re.match(r'.*k__(?P<kindom>\w+)(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?(s__(?P<species>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())
		# domain = m.groupdict()['domain']
		kindom = m.groupdict()['kindom']

		for sam in title:
			seqs += linedata[samindex[sam]]

			# samdomain[sam][domain] += linedata[samindex[sam]]
			samkindom[sam][kindom] += linedata[samindex[sam]]

			samphylum[sam][tmp[0]] += linedata[samindex[sam]]

			samclas[sam][tmp[1]] += linedata[samindex[sam]]

			samorder[sam][tmp[2]] += linedata[samindex[sam]]

			samfamily[sam][tmp[3]] += linedata[samindex[sam]]

			samgenus[sam][tmp[4]] += linedata[samindex[sam]]

			samall[sam] += linedata[samindex[sam]]

			samphylumper[sam]=defaultdict(float)
			samorderper[sam]=defaultdict(float)
			samfamilyper[sam]=defaultdict(float)
			samclasper[sam]=defaultdict(float)
			samgenusper[sam]=defaultdict(float)


		if tmp[0]=='Unclassified':
			for sam in title:
				samtaxstop[sam][domain] += linedata[samindex[sam]]
			

		if tmp[0] != 'Unclassified' and tmp[1] =='Unclassified':
			for sam in title:
				samtaxstop[sam][tmp[0]] += linedata[samindex[sam]]

		if tmp[1] != 'Unclassified' and tmp[2] =='Unclassified':
			for sam in title:
				samtaxstop[sam][tmp[1]] += linedata[samindex[sam]]

		if tmp[2] != 'Unclassified' and tmp[3] =='Unclassified':
			for sam in title:
				samtaxstop[sam][tmp[2]] += linedata[samindex[sam]]
		
		if tmp[3] != 'Unclassified' and tmp[4] =='Unclassified':
			for sam in title:
				samtaxstop[sam][tmp[3]] += linedata[samindex[sam]]
		

for sam in title:
	for sp in samphylum[sam].keys():
		samphylumper[sam][sp] = samphylum[sam][sp] / samall[sam]

	for sp in samclas[sam].keys():
		samclasper[sam][sp] = samclas[sam][sp] / samall[sam]

	for sp in samorder[sam].keys():
		samorderper[sam][sp] = samorder[sam][sp] / samall[sam]

	for sp in samfamily[sam].keys():
		samfamilyper[sam][sp] = samfamily[sam][sp] / samall[sam]

	for sp in samgenus[sam].keys():
		samgenusper[sam][sp] = samgenus[sam][sp] / samall[sam]

#p2c = defaultdict(list)
c2p = {}
#c2o = defaultdict(list)
o2c = {}
#o2f = defaultdict(list)
f2o = {}
#f2g = defaultdict(list)
g2f = {}

for p in taxline.keys():
#	p2c[p] = taxline[p].keys()
	for c in taxline[p].keys():
#		c2o[c] = taxline[p][c].keys()
		c2p[c] = p
		for o in taxline[p][c].keys():
#			o2f[o] = taxline[p][c][o].keys()
			o2c[o] = c
			for f in taxline[p][c][o].keys():
#				f2g[f] = taxline[p][c][o][f].keys()
				f2o[f] = o
				for g in taxline[p][c][o][f].keys():
					g2f[g] = f

path = options.out.rstrip('/')
for sam in title:

	out = open('{}/{}.tree.svg'.format(path, sam),'w')

	pickgenus=defaultdict(int)
	for sp in samgenus[sam].keys():	
		if not re.match(r'.*norank.*|.*uncultured.*|.*nclassified.*|.*Incertae_Sedis',sp):
			pickgenus[sp] += samgenus[sam][sp]

	pickgenussum = 0
	for i in pickgenus.keys():
		pickgenussum += pickgenus[i]

	pickgenusper = defaultdict(float)
	for i in pickgenus.keys():
		pickgenusper[i] = pickgenus[i]/pickgenussum


	gin=[]
	pick_genus = []
	
	hor = 300
	hormove = 240
	ver = 180
	vermove = 40
	r = 200
	samhi={}
	if options.config:
		config = open(options.config)
		for i in config:
			tmp = re.split(r'\s+',i.strip())
			if tmp[0]  not in title:
				print ('%s is not the exsit sample name' % tmp[0])
				sys.exit()
			samhi[tmp[0]] = int(tmp[1])
	if options.percent and not options.top:
		if sam not in samhi.keys():
			samhi[sam] = 3000
	
		print ('''<?xml version="1.0" encoding="utf-8"?>
	<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
	<svg version="1.1" width="2000" height="%s"  
	xmlns="http://www.w3.org/2000/svg"  
	xmlns:xlink="http://www.w3.org/1999/xlink"   
	xmlns:ev="http://www.w3.org/2001/xml-events"     
	baseProfile="full">''' % (samhi[sam]),file=out)
		for i in pickgenus.keys():
			if pickgenusper[i] >= options.percent:
				pick_genus.append(i)
				gin.append(i)
	elif options.top and not options.percent:
		print ('''<?xml version="1.0" encoding="utf-8"?>
	<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
	<svg version="1.1" width="2000" height="%s"  
	xmlns="http://www.w3.org/2000/svg"  
	xmlns:xlink="http://www.w3.org/1999/xlink"   
	xmlns:ev="http://www.w3.org/2001/xml-events"     
	baseProfile="full">''' % (2000 + ver*(options.top-10)),file=out)
		sorted_genus = sorted(pickgenus.items(),key=lambda x:x[1],reverse=True)
		if options.top <= len(sorted_genus):
			pick_genus = sorted_genus[0:options.top]
		else:
			pick_genus = sorted_genus
		for i in pick_genus:
			gin.append(i[0])

	else:
		print ('Warning:  you must specify -t or -p, if choose -p, you must specify the "--hi" to make a proper height ')
		sys.exit()
	fin = defaultdict(list)
	for i in gin:
		fin[g2f[i]].append(i)

	oin = defaultdict(list)
	for i in fin.keys():
		oin[f2o[i]].append(i)

	cin = defaultdict(list)
	for i in oin.keys():
		cin[o2c[i]].append(i)

	pin = defaultdict(list)
	for i in cin.keys():
		pin[c2p[i]].append(i)

	spper = {}
	
	outfile = open('{}/{}.tax_percent.txt'.format(path, sam),'w')
	print ('Tax\t%s' % sam,file = outfile)
	print ('phylum',file = outfile)
	for sp in pin.keys():
		print (sp,samphylum[sam][sp],'%.3f' % samphylumper[sam][sp],sep = '\t',file = outfile)
	print ('\nclass',file = outfile)
	for sp in cin.keys():
		print (sp,samclas[sam][sp],'%.3f' % samclasper[sam][sp],sep='\t',file = outfile)
	print ('\norder',file = outfile)
	for sp in oin.keys():
		print (sp,samorder[sam][sp],'%.3f' % samorderper[sam][sp],sep='\t',file = outfile)
	print ('\nfamily',file = outfile)
	for sp in fin.keys():
		print (sp,samfamily[sam][sp],'%.3f' % samfamilyper[sam][sp],sep='\t',file=outfile)
	print ('\ngenus',file=outfile)
	for sp in gin:
		print (sp,samgenus[sam][sp],'%.3f' % samgenusper[sam][sp],sep='\t',file=outfile)




	dcoord = defaultdict(list)
	pcoord = defaultdict(list)
	ccoord = defaultdict(list)
	ocoord = defaultdict(list)
	fcoord = defaultdict(list)
	gcoord = defaultdict(list)
	n = 2
	py=[]
	for p in sorted(pin.keys()):
		cy = []
		for c in sorted(pin[p]):
			oy = []
			for o in sorted(cin[c]):
				fy = []
				for f in sorted(oin[o]):
					gy = []
					for g in sorted(fin[f]):
						gcoord[g] = [5*hor+hormove, ver*(n-1)+vermove]
						n += 1
						gy.append(gcoord[g][1])
	
					fcoord[f] = [4*hor+hormove, (max(gy)+min(gy))/2]
					fy.append(fcoord[f][1])
	
				ocoord[o] = [3*hor+hormove, (max(fy)+min(fy))/2]
				oy.append(ocoord[o][1])
	
			ccoord[c] = [2*hor+hormove, (max(oy)+min(oy))/2]
			cy.append(ccoord[c][1])
	
		pcoord[p] = [1*hor+hormove, (max(cy)+min(cy))/2]
		py.append(pcoord[p][1])

	try:
		dcoord=[hormove,(max(py)+min(py))/2]
	except:
		sys.stderr.write('\nError: no species abundance percentage satisfiy the cutoff!\n\n')
		sys.exit(1)

########改##########
	## 添加颜色
	color=['#EE3B3B','#458B00','#FFD700','#8B7D6B','#66CD00','#0000FF','#8A2BE2','#B8860B','#FF82AB','#6495ED']
	tax_list = ['phylum', 'class', 'order', 'family', 'genus']
	lhi = 0
	vermove = 100
	for i, j in enumerate(tax_list):
		print ('<rect x="80" y="%s" width="70" height="20" style="fill:%s"/>' % (vermove+lhi+100,color[i]),file=out)
		print ('<text x="165" y="%s" fill="#000000" font-size="20" alignment-baseline="hanging">%s</text>' % (vermove+lhi+100+14,j),file=out)
		lhi += 35
########改##########


#print (dcoord[domain][0],dcoord[domain][1])
	print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (dcoord[0],dcoord[1],dcoord[0]+hor/2,dcoord[1]),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (dcoord[0],dcoord[1]-25,domain),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (dcoord[0],dcoord[1],3),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (dcoord[0],dcoord[1]+25,samtaxstop[sam][i],seqs),file=out)


	for p in sorted(pin.keys()):
		print ('<polyline points="%s,%s %s,%s %s,%s" style="fill:none;stroke:black;stroke-width:2"/>' % (dcoord[0]+hor/2,dcoord[1],dcoord[0]+hor/2,pcoord[p][1],pcoord[p][0],pcoord[p][1]),file=out)
		if len(pin[p]) >1:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (pcoord[p][0],pcoord[p][1],  pcoord[p][0]+hor/2,pcoord[p][1]) ,file =out)
			for c in sorted(pin[p]):
				print ('<polyline points="%s,%s %s,%s %s,%s" style="fill:none;stroke:black;stroke-width:2"/>' % (pcoord[p][0]+hor/2,pcoord[p][1],pcoord[p][0]+hor/2,ccoord[c][1],ccoord[c][0],ccoord[c][1]),file=out)
		else:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (pcoord[p][0],pcoord[p][1],ccoord[pin[p][0]][0],ccoord[pin[p][0]][1]),file=out)



	for c in sorted(cin.keys()):
		if len(cin[c]) >1:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (ccoord[c][0],ccoord[c][1],  ccoord[c][0]+hor/2,ccoord[c][1]) ,file =out)
			for o in cin[c]:
				print ('<polyline points="%s,%s %s,%s %s,%s" style="fill:none;stroke:black;stroke-width:2"/>' % (ccoord[c][0]+hor/2,ccoord[c][1],ccoord[c][0]+hor/2,ocoord[o][1],ocoord[o][0],ocoord[o][1]),file=out)
		else:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (ccoord[c][0],ccoord[c][1],ocoord[cin[c][0]][0],ocoord[cin[c][0]][1]),file=out)
		
			

	for o in oin.keys():
		if len(oin[o]) >1:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (ocoord[o][0],ocoord[o][1],  ocoord[o][0]+hor/2,ocoord[o][1]) ,file =out)
			for f in oin[o]:
				print ('<polyline points="%s,%s %s,%s %s,%s" style="fill:none;stroke:black;stroke-width:2"/>' % (ocoord[o][0]+hor/2,ocoord[o][1],ocoord[o][0]+hor/2,fcoord[f][1],fcoord[f][0],fcoord[f][1]),file=out)
		else:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (ocoord[o][0],ocoord[o][1],fcoord[oin[o][0]][0],fcoord[oin[o][0]][1]),file=out)


	for f in sorted(fin.keys()):
		if len(fin[f]) >1:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (fcoord[f][0],fcoord[f][1],  fcoord[f][0]+hor/2,fcoord[f][1]) ,file =out)
			for g in fin[f]:
				print ('<polyline points="%s,%s %s,%s %s,%s" style="fill:none;stroke:black;stroke-width:2"/>' % (fcoord[f][0]+hor/2,fcoord[f][1],fcoord[f][0]+hor/2,gcoord[g][1],gcoord[g][0],gcoord[g][1]),file=out)
		else:
			print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (fcoord[f][0],fcoord[f][1],gcoord[fin[f][0]][0],gcoord[fin[f][0]][1]),file=out)



	for i in pcoord.keys():
		R = samphylumper[sam][i]*r
		if R < 20:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (pcoord[i][0],pcoord[i][1]-25,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (pcoord[i][0],pcoord[i][1]+28,samtaxstop[sam][i],samphylum[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#EE3B3B"/>' % (pcoord[i][0],pcoord[i][1],samphylumper[sam][i]*r),file=out)
		else:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (pcoord[i][0],pcoord[i][1]-R-5,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (pcoord[i][0],pcoord[i][1]+R+20,samtaxstop[sam][i],samphylum[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#EE3B3B"/>' % (pcoord[i][0],pcoord[i][1],samphylumper[sam][i]*r),file=out)
		

	for i in ccoord.keys():
		R = samclasper[sam][i]*r
		if R < 20:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ccoord[i][0],ccoord[i][1]-25,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ccoord[i][0],ccoord[i][1]+30,samtaxstop[sam][i],samclas[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#458B00"/>' % (ccoord[i][0],ccoord[i][1],samclasper[sam][i]*r),file=out)
		else:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ccoord[i][0],ccoord[i][1]-R-5,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ccoord[i][0],ccoord[i][1]+R+20,samtaxstop[sam][i],samclas[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#458B00"/>' % (ccoord[i][0],ccoord[i][1],samclasper[sam][i]*r),file=out)
		

	for i in ocoord.keys():
		R = samorderper[sam][i]*r
		if R < 20:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ocoord[i][0],ocoord[i][1]-25,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ocoord[i][0],ocoord[i][1]+30,samtaxstop[sam][i],samorder[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#FFD700"/>' % (ocoord[i][0],ocoord[i][1],samorderper[sam][i]*r),file=out)
		else:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ocoord[i][0],ocoord[i][1]-5-R,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ocoord[i][0],ocoord[i][1]+20+R,samtaxstop[sam][i],samorder[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#FFD700"/>' % (ocoord[i][0],ocoord[i][1],samorderper[sam][i]*r),file=out)
		

	for i in fcoord.keys():
		R = samfamilyper[sam][i]*r
		if R < 20:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (fcoord[i][0],fcoord[i][1]-25,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (fcoord[i][0],fcoord[i][1]+30,samtaxstop[sam][i],samfamily[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#8B7D6B"/>' % (fcoord[i][0],fcoord[i][1],samfamilyper[sam][i]*r),file=out)
		else:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (fcoord[i][0],fcoord[i][1]-5-R,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (fcoord[i][0],fcoord[i][1]+20+R,samtaxstop[sam][i],samfamily[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#8B7D6B"/>' % (fcoord[i][0],fcoord[i][1],samfamilyper[sam][i]*r),file=out)
		

	for i in gcoord.keys():
		R = samgenusper[sam][i]*r
		if R < 20:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (gcoord[i][0],gcoord[i][1]-25,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (gcoord[i][0],gcoord[i][1]+30,samtaxstop[sam][i],samgenus[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#66CD00"/>' % (gcoord[i][0],gcoord[i][1],samgenusper[sam][i]*r),file=out)
		else:
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (gcoord[i][0],gcoord[i][1]-5-R,i),file=out)
			print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (gcoord[i][0],gcoord[i][1]+20+R,samtaxstop[sam][i],samgenus[sam][i]),file=out)
			print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="#66CD00"/>' % (gcoord[i][0],gcoord[i][1],samgenusper[sam][i]*r),file=out)
		

	print ('</svg>',file=out)

	out.close()
	outfile.close()

	os.system('convert {}/{}.tree.svg {}/{}.tree.pdf'.format(path, sam, path, sam))
	os.system('convert {}/{}.tree.svg {}/{}.tree.png'.format(path, sam, path, sam))













