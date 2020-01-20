from collections import defaultdict
import re
#d=tree()
#file = open('otu_taxa_table.xls')
#line=[]
def levelline(i):
		def tree(): return defaultdict(tree)
		d = tree()
		# m = re.match(r'.*d__(?P<domain>\w+)(;\s)?(k__(?P<kindom>((\w+|-|\.)+)))?(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())
		m = re.match(r'.*k__(?P<kindom>\w+)(;\s)?(p__(?P<phylum>((\w+|-|\.)+)))?(;\s)?(c__(?P<class>((\w|-|\.)+)))?(;\s)?(o__(?P<order>((\w|-|\.)+)))?(;\s)?(f__(?P<family>((\w|-|\.)+)))?(;\s)?(g__(?P<genus>((\w|-|\.)+)))?(;\s)?(s__(?P<species>((\w|-|\.)+)))?(;\s)?.*' ,i.strip())

		#print m
		#print i
		if m:
			if not m.groupdict()['phylum']:
				p = 'Unclassified'
			else:
				if re.match(r'norank',m.groupdict()['phylum']):
					if not m.groupdict()['kindom']:
						p = m.groupdict()['domain'] + '_norank'
					else:
						p = m.groupdict()['kindom'] + '_norank'
				else:
					p = m.groupdict()['phylum']
	
	
			if not m.groupdict()['class']:
	#			m.groupdict()['class']='Unclassified'
				c = p + '_unclassified'
			else:
				if re.match(r'norank',m.groupdict()['class']):
					if not re.match(r'.*norank.*',p):
						c = p + '_norank'
					else:
						c = p
				elif re.match(r'.*Incertae_Sedis',m.groupdict()['class']):
					tmp = re.match(r'(Class)?(.*Incertae_Sedis)',m.groupdict()['class'])
					if not re.match(r'.*Incertae_Sedis',p):
						c = p + '_' + tmp.group(2).lower()
				#	else:
				#		c = p
				elif re.match(r'Unknown_Class',m.groupdict()['class']):
					c = p + '_unknown'
				else:
					c = m.groupdict()['class']
	
			if not m.groupdict()['order']:
				o = 'Unclassified'
			else:
				if re.match(r'norank',m.groupdict()['order']):
					if not re.match(r'.*norank.*',c):
						o = c +'_norank'
					else:
						o = c
				elif re.match(r'.*Incertae_Sedis',m.groupdict()['order']):
					tmp = re.match(r'(Order)?(.*Incertae_Sedis)',m.groupdict()['order'])
					if not re.match(r'.*Incertae_Sedis',c):
						o = c + tmp.group(2)
					elif not re.match(r'.*Incertae_Sedis',p):
						o = p + tmp.group(2)
				#	else:
				#		o = c
				elif re.match(r'Unknown_Order',m.groupdict()['order']):
					o = c + 'unknown'
				else:
					o = m.groupdict()['order']
	
			if not m.groupdict()['family']:
				f = 'Unclassified'
			else:
				if re.match(r'norank',m.groupdict()['family']):
					if not re.match(r'.*norank.*' , o):
						f = o + '_norank'
					else:
						f = o
				elif re.match(r'Family_.*Incertae_Sedis',m.groupdict()['family']):
					tmp = re.match(r'(Family_.*Incertae_Sedis)',m.groupdict()['family'])
#					if not re.match(r'.*Incertae_Sedis',o):
					f = o + "_" + tmp.group(1)
#					elif not re.match(r'.*Incertae_Sedis',c):
#						f = c + tmp.group(1)
#					elif not re.match(r'.*Incertae_Sedis',p):
#						f = p + tmp.group(1)
	#				else:
	#					f = o
				elif re.match(r'FamilyI+',m.groupdict()['family']):
					tmp = re.match(r'(Family(I|V|X)+)',m.groupdict()['family'])
					f = o + '_' + tmp.group(1)
				elif re.match(r'Incertae_Sedis' , m.groupdict()['family']):
					f = o + '_' + 'Incertae_Sedis'
				elif re.match(r'Unknown_Family',m.groupdict()['family']):
					f = o + 'unknown'
				else:
					f = m.groupdict()['family']
	
			if not m.groupdict()['genus']:
				g = 'Unclassified'
			else:
				if re.match(r'norank',m.groupdict()['genus']):
					if not re.match(r'.*norank.*' , f):
						g = f + '_norank'
					else:
						g = f
				elif re.match(r'.*Incertae_Sedis',m.groupdict()['genus']):
					tmp = re.match(r'(Genus)?(.*Incertae_Sedis)',m.groupdict()['genus'])
					if not re.match(r'.*Incertae_Sedis',f):
						g = f + '_' + tmp.group(2)
					elif not re.match(r'.*Incertae_Sedis',o):
						g = o + '_' + tmp.group(2)
					elif not re.match(r'.*Incertae_Sedis',c):
						g = c + '_' + tmp.group(2)
					elif not re.match(r'.*Incertae_Sedis',p):
						g = p + '_' + tmp.group(2)
				#	else:
				#		g = f
				elif re.match(r'Unknown_Genus',m.groupdict()['genus']):
					g = f + 'unknown'
				else:
					g = m.groupdict()['genus']

			if d[p][c][o][f][g]==1:
				pass
			else:
				d[p][c][o][f][g] = 1
				return (p,c,o,f,g)
			
