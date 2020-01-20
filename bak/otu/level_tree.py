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
p.add_option('-g',dest='group',help='specify how to fen zu')
p.add_option('-o',dest='out',help='the path of out')
p.add_option('--hi',dest='height',help="specify the height of figure, suggest 2000+ if species > 10",type="int")
(options,args) = p.parse_args()


def tree(): return defaultdict(tree)

file = open(options.file)
taxline = tree()
phylum = defaultdict(int)
clas =  defaultdict(int)
order = defaultdict(int)
family = defaultdict(int)
genus = defaultdict(int)

taxstop = defaultdict(int)

samindex = defaultdict(int)

samindex = {}
samphylum =	{}
samclas = {}
samorder = {}
samfamily = {}
samgenus = {}
samkindom = {}

seqs = 0
fortop=[]
title = re.split(r'\t',file.readline().strip())[1:-1]

gp = defaultdict(list)
picksample=[]
if options.group:
	groupfile = open(options.group)
	for i in groupfile:
		tmp = re.split(r'\s+', i.strip())
		gp[tmp[1]].append(tmp[0])
		picksample.append(tmp[0])
else:
	for i in title:
		gp[i].append(i)
	picksample = title


for i in title:
	samindex[i] = title.index(i)
	samphylum[i] = defaultdict(int)
	samclas[i] = defaultdict(int)
	samorder[i] = defaultdict(int)
	samfamily[i] = defaultdict(int)
	samgenus[i] = defaultdict(int)
	samkindom[i] = defaultdict(int)

samphylumall=defaultdict(int)
samclasall=defaultdict(int)
samorderall=defaultdict(int)
samfamilyall=defaultdict(int)
samgenusall=defaultdict(int)
kindom=''
for i in file.readlines():
	tmp = levelline(i.strip())
	#print(tmp)
	if tmp:
		taxline[tmp[0]][tmp[1]][tmp[2]][tmp[3]][tmp[4]] = 1
		linedata = [int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]]

		for sam in picksample:
			seqs += linedata[samindex[sam]]

			samkindom[sam][kindom] += linedata[samindex[sam]]

			samphylum[sam][tmp[0]] += linedata[samindex[sam]]
			samphylumall[tmp[0]] += linedata[samindex[sam]]

			samclas[sam][tmp[1]] += linedata[samindex[sam]]
			samclasall[tmp[1]] += linedata[samindex[sam]]

			samorder[sam][tmp[2]] += linedata[samindex[sam]]
			samorderall[tmp[2]] += linedata[samindex[sam]]

			samfamily[sam][tmp[3]] += linedata[samindex[sam]]
			samfamilyall[tmp[3]] += linedata[samindex[sam]]

			samgenus[sam][tmp[4]] += linedata[samindex[sam]]
			samgenusall[tmp[4]] += linedata[samindex[sam]]
		
		# m = re.match(r'.*(d__(?P<domain>\w+))(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())
		m = re.match(r'.*k__(?P<kindom>\w+)(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?(s__(?P<species>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())
		kindom = m.groupdict()['kindom']

		if tmp[0]=='Unclassified':
			for j in gp.keys():
				for k in gp[j]:
					taxstop[kindom] += int(float(re.split(r'\t',i.strip())[1:-1][samindex[k]]))

		if tmp[0] != 'Unclassified' and tmp[1] =='Unclassified':
			for j in gp.keys():
				for k in gp[j]:
					taxstop[tmp[0]] += int(float(re.split(r'\t',i.strip())[1:-1][samindex[k]]))

		if tmp[1] != 'Unclassified' and tmp[2] =='Unclassified':
			for j in gp.keys():
				for k in gp[j]:
					taxstop[tmp[1]] += int(float(re.split(r'\t',i.strip())[1:-1][samindex[k]]))

		if tmp[2] != 'Unclassified' and tmp[3] =='Unclassified':
			for j in gp.keys():
				for k in gp[j]:
					taxstop[tmp[2]] += int(float(re.split(r'\t',i.strip())[1:-1][samindex[k]]))
		
		if tmp[3] != 'Unclassified' and tmp[4] =='Unclassified':
			for j in gp.keys():
				for k in gp[j]:
					taxstop[tmp[3]] += int(float(re.split(r'\t',i.strip())[1:-1][samindex[k]]))
			

		phylum[tmp[0]] += sum([int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]])
		clas[tmp[1]] += sum([int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]])
		order[tmp[2]] += sum([int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]])
		family[tmp[3]] += sum([int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]])
		if not re.match(r'.*norank.*|.*uncultured.*|.*Unclassified.*|.*Incertae_Sedis',tmp[4]):
			genus[tmp[4]] += sum([int(float(x)) for x in re.split(r'\t',i.strip())[1:-1]])


gpkindom = {}
gpphylum = {}
gpclas = {}
gporder = {}
gpfamily = {}
gpgenus = {}


for i in gp.keys():
	gpphylum[i] = defaultdict(int)
	gpclas[i] = defaultdict(int)
	gporder[i] = defaultdict(int)
	gpfamily[i] = defaultdict(int)
	gpgenus[i] = defaultdict(int)
	gpkindom[i] = defaultdict(int)

for grp in gp.keys():
	for sp in phylum.keys():
		for sam in gp[grp]:
			gpphylum[grp][sp] += samphylum[sam][sp]

	for sp in clas.keys():
		for sam in gp[grp]:
			gpclas[grp][sp] += samclas[sam][sp]

	for sp in order.keys():
		for sam in gp[grp]:
			gporder[grp][sp] += samorder[sam][sp]
	
	for sp in family.keys():
		for sam in gp[grp]:
			gpfamily[grp][sp] += samfamily[sam][sp]

	for sp in genus.keys():
		for sam in gp[grp]:
			gpgenus[grp][sp] += samgenus[sam][sp]
	
	for sam in gp[grp]:
		gpkindom[grp][kindom] += samkindom[sam][kindom]

pickgenus=defaultdict(int)
for sam in picksample:
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
path = options.out.rstrip('/')
out = open('{}/level_tree.svg'.format(path),'w')

hor = 300
hormove = 240
ver = 180
vermove = 40


if options.percent and not options.top:
	

	print ('''<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg version="1.1" width="2000" height="%s"  
xmlns="http://www.w3.org/2000/svg"  
xmlns:xlink="http://www.w3.org/1999/xlink"   
xmlns:ev="http://www.w3.org/2001/xml-events"     
baseProfile="full">''' % (options.height),file=out)
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
	pick_genus = sorted_genus[0:options.top]
	for i in pick_genus:
		gin.append(i[0])
else:
	print ('Warning:  you must specify -t or -p, if choose -p, you must specify the "--hi" to make a proper height ')
	sys.exit(1)


p2c = defaultdict(list)
c2p = {}
c2o = defaultdict(list)
o2c = {}
o2f = defaultdict(list)
f2o = {}
f2g = defaultdict(list)
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

outfile = open('{}/tax_percent.xls'.format(path),'w')
print ('Tax\t%s' % '\t'.join(sorted(gp.keys())),file = outfile)

spper = {}


print (kindom,end='\t',file=outfile)
spper[kindom] = defaultdict(int)
dosum = 0
for grp in gp.keys():
	dosum += gpkindom[grp][kindom]
for grp in sorted(gp.keys()):
	spper[kindom][grp] = gpkindom[grp][kindom]/dosum
	print ('%.3f' % spper[kindom][grp],end='\t',file=outfile)
print ('\n',file=outfile)

for sp in pin.keys():
	print (sp,end='\t',file=outfile)
	spper[sp] = defaultdict(list)
	spsum = 0
	for grp in gp.keys():
		spsum += gpphylum[grp][sp]
	for grp in sorted(gp.keys()):
		spper[sp][grp] = gpphylum[grp][sp]/spsum
		print ('%.3f' % spper[sp][grp],end='\t',file=outfile)
	print ('',file=outfile)
print ('',file=outfile)
for sp in cin.keys():
	print (sp,end='\t',file=outfile)
	spper[sp] = defaultdict(list)
	spsum = 0
	for grp in gp.keys():
		spsum += gpclas[grp][sp]
	for grp in gp.keys():
		spper[sp][grp] = gpclas[grp][sp]/spsum
		print ('%.3f' % spper[sp][grp],end='\t',file=outfile)
	print ('',file=outfile)
#
print ('',file=outfile)
for sp in oin.keys():
	print (sp,end='\t',file=outfile)
	spper[sp] = defaultdict(list)
	spsum = 0
	for grp in gp.keys():
		spsum += gporder[grp][sp]
	for grp in gp.keys():
		spper[sp][grp] = gporder[grp][sp]/spsum
		print ('%.3f' % spper[sp][grp],end='\t',file=outfile)
	print ('',file=outfile)

print ('',file=outfile)
for sp in fin.keys():
	print (sp,end='\t',file=outfile)
	spper[sp] = defaultdict(list)
	spsum = 0
	for grp in gp.keys():
		spsum += gpfamily[grp][sp]
	for grp in gp.keys():
		spper[sp][grp] = gpfamily[grp][sp]/spsum
		print ('%.3f' % spper[sp][grp],end='\t',file=outfile)
	print ('',file=outfile)

print ('',file=outfile)
for sp in gin:
	print (sp,end='\t',file=outfile)
	spper[sp] = defaultdict(list)
	spsum = 0
	for grp in gp.keys():
		spsum += gpgenus[grp][sp]
	for grp in gp.keys():
		spper[sp][grp] = gpgenus[grp][sp]/spsum
		print ('%.3f' % spper[sp][grp],end='\t',file=outfile)
	print ('',file=outfile)


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

color=['#EE3B3B','#458B00','#FFD700','#8B7D6B','#66CD00','#0000FF','#8A2BE2','#B8860B','#FF82AB','#6495ED', '#00CC99', '#99CCCC']
gpcolor={}

try:
	for i in sorted(gp.keys()):
		gpcolor[i] = color.pop(0)
except:
	sys.stderr.write('\nError: Samples or Groups Number Must be equal or less then 12!\n\n')
	sys.exit(1)



lhi = 0
for i in sorted(gpcolor.keys()):
	print ('<rect x="80" y="%s" width="70" height="20" style="fill:%s"/>' % (vermove+lhi+100,gpcolor[i]),file=out)
	print ('<text x="165" y="%s" fill="#000000" font-size="20" alignment-baseline="hanging">%s</text>' % (vermove+lhi+100+14,i),file=out)
	lhi += 35

#print (dcoord[kindom][0],dcoord[kindom][1])
print ('<polyline points="%s,%s %s,%s " style="fill:none;stroke:black;stroke-width:2"/>' % (dcoord[0],dcoord[1],dcoord[0]+hor/2,dcoord[1]),file=out)
print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (dcoord[0],dcoord[1]+28,kindom),file=out)
print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (dcoord[0],dcoord[1],3),file=out)
print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (dcoord[0],dcoord[1]+56,taxstop[i],seqs),file=out)
for i in [kindom]:
	start = [dcoord[0]+40, dcoord[1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [dcoord[0]+40, dcoord[1]-60]
		else:
			stop = [dcoord[0]+cos(2*pi*spper[i][grp]+tmp)*40, dcoord[1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (dcoord[0],dcoord[1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
#		print (i,grp,start,stop)
		start = stop

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
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (pcoord[i][0],pcoord[i][1]+28,i),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (pcoord[i][0],pcoord[i][1]+56,taxstop[i],samphylumall[i]),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (pcoord[i][0],pcoord[i][1],3),file=out)
	start = [pcoord[i][0]+40, pcoord[i][1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
	#	print (i,grp,spper[i][grp])
		

		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [pcoord[i][0]+40, pcoord[i][1]-60]
		else:
			stop = [pcoord[i][0]+cos(2*pi*spper[i][grp]+tmp)*40, pcoord[i][1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (pcoord[i][0],pcoord[i][1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
#		print (i,grp,start,stop)
		start = stop

for i in ccoord.keys():
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ccoord[i][0],ccoord[i][1]+28,i),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ccoord[i][0],ccoord[i][1]+56,taxstop[i],samclasall[i]),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (ccoord[i][0],ccoord[i][1],3),file=out)
	start = [ccoord[i][0]+40, ccoord[i][1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [ccoord[i][0]+40, ccoord[i][1]-60]
		else:
			stop = [ccoord[i][0]+cos(2*pi*spper[i][grp]+tmp)*40, ccoord[i][1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (ccoord[i][0],ccoord[i][1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
		start = stop
	

for i in ocoord.keys():
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (ocoord[i][0],ocoord[i][1]+28,i),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (ocoord[i][0],ocoord[i][1]+56,taxstop[i],samorderall[i]),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (ocoord[i][0],ocoord[i][1],3),file=out)
	start = [ocoord[i][0]+40, ocoord[i][1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [ocoord[i][0]+40, ocoord[i][1]-60]
		else:
			stop = [ocoord[i][0]+cos(2*pi*spper[i][grp]+tmp)*40, ocoord[i][1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (ocoord[i][0],ocoord[i][1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
		start = stop

for i in fcoord.keys():
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (fcoord[i][0],fcoord[i][1]+28,i),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (fcoord[i][0],fcoord[i][1]+56,taxstop[i],samfamilyall[i]),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (fcoord[i][0],fcoord[i][1],3),file=out)
	start = [fcoord[i][0]+40, fcoord[i][1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [fcoord[i][0]+40, fcoord[i][1]-60]
		else:
			stop = [fcoord[i][0]+cos(2*pi*spper[i][grp]+tmp)*40, fcoord[i][1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (fcoord[i][0],fcoord[i][1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
		start = stop

for i in gcoord.keys():
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="25">%s</text>' % (gcoord[i][0],gcoord[i][1]+28,i),file=out)
	print ('<text x="%s" y="%s" fill="#000000" text-anchor="middle" font-size="22">%s , %s</text>' % (gcoord[i][0],gcoord[i][1]+56,taxstop[i],samgenusall[i]),file=out)
	print ('<circle cx="%s" cy="%s" r="%s" stroke="black" stroke-width="2" fill="black"/>' % (gcoord[i][0],gcoord[i][1],3),file=out)
	start = [gcoord[i][0]+40, gcoord[i][1]-60]
	tmp = 0
	for grp in sorted(gp.keys()):
		if sorted(gp.keys()).index(grp) == len(gp.keys())-1:
			stop = [gcoord[i][0]+40, gcoord[i][1]-60]
		else:
			stop = [gcoord[i][0]+cos(2*pi*spper[i][grp]+tmp)*40, gcoord[i][1]-60-sin(2*pi*spper[i][grp]+tmp)*40]
			tmp += 2*pi*spper[i][grp]
		if 2*pi*spper[i][grp] < pi:
			hu = 0
		else:
			hu = 1
		print ('<path d="M%s,%s  L%s,%s  A%s,%s 0 %s,0 %s,%s z" style="fill:%s"/>' % (gcoord[i][0],gcoord[i][1]-60 ,start[0],start[1], 40,40, hu ,stop[0],stop[1], gpcolor[grp]),file=out)
		start = stop

print ('</svg>',file=out)

out.close()

os.system('convert {}/level_tree.svg {}/level_tree.png'.format(path, path))
os.system('convert {}/level_tree.svg {}/level_tree.pdf'.format(path, path))














