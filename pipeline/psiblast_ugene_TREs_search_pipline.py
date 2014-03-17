import os
import sys
from Bio.Blast.Applications import NcbitblastnCommandline

                                                                                                        ## denotes changable run setting
#establish work directory
path_to_work = sys.argv[5]                                                                              ##
work_directory = sys.argv[1]                                                                            ##
work_path = path_to_work + work_directory + '/'
print 'Results in ' + work_path
if not os.path.exists(work_path):
        os.makedirs(work_path)
else:
    raise Exception('work directory already exists')
os.makedirs(work_path + 'fragments_diagrams')
os.makedirs(work_path + 'TEs_diagrams')
pssm_directory = sys.argv[6]

#establish genome
genome = sys.argv[2] + sys.argv[3]                                                                      ##
genome_code = work_directory                                                                            ##

#make blast db
print 'Assembly is in ' + genome
os.system('makeblastdb -in ' + genome + ' -dbtype nucl')
outnamestart = work_path + genome_code + '_' #This is the work directory + the genome code
                                             #to establish the start of output file names
                                             #(e.g. /work_path/Ppac_)

#establish YR fragment file. The ending will be added depanding on the output format (e.g. fasta, embl, etc..)
fragments_file = outnamestart + 'YR_fragments'

#Other parameters
extend = 10000    #The number of BP to extend the YR hit in both directions                             ##
intron = 100 #The length of sequence between hits of the same domain that can be considered as intron   ##
print 'extend = ' + str(extend)
print 'max intron = ' + str(intron)

#blast parameters
e_value = 0.01    #For keeping blast results                                                            ##
print 'e value cutoff = ' + str(e_value)

#ugene parameters
ident = 100       #minimum identity between repeat partners                                             ##
maxdist = 20000   #maximum distance between repeat partners                                             ##
mindist = 0       #minimum distance between repeat partners                                             ##
minlength = 20    #minimum repeat length                                                                ##
comment = sys.argv[4]                  ##
                                                                                                        ##No changable settings after here
print 'Eugene parameters: identity=' + str(ident) + ' maxdist=' + str(maxdist) + ' mindist=' + str(mindist) + ' minlength=' + str(minlength)
print comment




#Run psitblastn  ######################################################
def run_tblastn(pssm_file, d_b, e_value, outnamestart):
    pssm_name = pssm_file.split('/')[-1].split('.')[0] # pssm name is either YR, RT or MT
    tblastn_cline = NcbitblastnCommandline(in_pssm = pssm_file, db=d_b, evalue=e_value, outfmt=5, out=outnamestart + pssm_name + '.xml')
    stdout, stderr = tblastn_cline()
    return outnamestart + pssm_name + '.xml'
#######################################################################




#Run YR psitblastn on the genome assembly  ############################
#print 'Searching YR on genome'
pssm_file = pssm_directory + 'YR.ckp'
d_b = genome
YR_xml = run_tblastn(pssm_file, d_b, e_value, outnamestart)
#######################################################################




#Make YR fragments fasta and YR fragments SeqRecord object. 
#YR fragments are YR hits extanded by the number of BPs that are in
#'extend' in both dirrections

from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
import re
import sys




#Get a contig sequence given contig name  ##############################
recs = SeqIO.index(genome, "fasta")
def seq_for_contig(contig_name,recs):
    return recs[contig_name].seq
########################################################################



YR_blast_qresult = SearchIO.read(YR_xml, 'blast-xml')
if len(YR_blast_qresult) == 0:
    raise Exception('No YR were found')

#make fragments
#print 'Preparing YR fragments'
fragment_stores = []
for hit in YR_blast_qresult:
    hit_sequence = seq_for_contig(hit.id, recs)
    for hsp in hit.hsps:
        fragstart = 0
	fragend = len(hit_sequence)
	if hsp.hit_start > extend:
	    fragstart = hsp.hit_start - extend
	if fragend - hsp.hit_end > extend:
	    fragend =  hsp.hit_end + extend
        hsp_id = hit.id + ':' + str(fragstart) + ':' + str(fragend)
        hsp_sequence = Seq(str(hit_sequence[fragstart:fragend]),IUPAC.ambiguous_dna)
	fragment_record = SeqRecord(hsp_sequence, id = hsp_id + ":YR_start_" + str(hsp.hit_start))
            #print '  Got fragment ' + fragment_record.id
	fragment_stores.append(fragment_record)

SeqIO.write(fragment_stores, fragments_file + '.fas', "fasta")

#Turn fragments file to a blast DB
#print 'Making fragments blast db'

os.system('makeblastdb -in ' + fragments_file + '.fas -dbtype nucl')




#Run RT psitblastn on the YR fragments#################################
#print 'Searching RT in fragments'
pssm_file = pssm_directory + 'RT.ckp'
d_b = fragments_file + '.fas'
RT_xml = run_tblastn(pssm_file, d_b, e_value, outnamestart)
#######################################################################





#Run MT psitblastn on the YR fragments#################################
#print 'Searching MT in fragments'
pssm_file = pssm_directory + 'MT.ckp'
MT_xml = run_tblastn(pssm_file, d_b, e_value, outnamestart)
#######################################################################




RT_blast_qresult = SearchIO.read(RT_xml, 'blast-xml')
MT_blast_qresult = SearchIO.read(MT_xml, 'blast-xml')
 




#Add domain annotations to fragment sequence





#Make DOMAIN feature#######################################################
def make_feature(product, blast_qresult, fragment ,hit, hsp, fragstart, count):
    s = hsp.hit_start
    e = hsp.hit_end
    if product == 'YR':
        s = hsp.hit_start-int(fragstart)
        e = hsp.hit_end-int(fragstart)
    feature = SeqFeature(FeatureLocation(s, e), type="DOMAIN", strand= hsp.hit_strand)
    feature.qualifiers['loc_on_contig'] = str(hsp.hit_start+1) + '..' + str(hsp.hit_end)
    feature.qualifiers['product'] = product
    feature.qualifiers['serial_on_frag'] = count
    count += 1
    feature.qualifiers['program'] = blast_qresult.program + "_" + blast_qresult.version
    feature.qualifiers['evalue'] = hsp.evalue
    feature.qualifiers['assembly'] = blast_qresult.target.split('/')[-1]
    feature.qualifiers['contig'] = contig
    feature.qualifiers['translation'] = feature.extract(fragment.seq).translate()
    return (feature, count)
##########################################################################




#print 'Annotating features to fragments'                    
for fragment in fragment_stores:
    
    #Modify the fragment names
    keep_original_id = fragment.id
    contig = fragment.id.split(':')[0]
    fragstart = fragment.id.split(':')[1]
    fragend = fragment.id.split(':')[2]
    YR_start = fragment.id.split(':')[3]
    hsp_identifier = fragment.id.split('_')[-1]
    #Correct fragment description
    fragment.description = 'Assembly: ' + YR_blast_qresult.target.split('/')[-1] + ', contig: ' + contig + ' from: ' + str(int(fragstart)+1) + ' to: ' + fragend
    #Correct fragment ID
    identification = contig
    for i in ('_','-',';',':','/'):
        if i in contig:
            identification = contig.split(i)[-1]
    ID = re.sub(r'\D','',identification) 
    identification = ID + '_' + str(int(hsp_identifier) + 1)
    fragment.id = identification
    
    #Set fragment annotations
    fragment.annotations['fragment_starts_at'] = int(fragstart)+1
    fragment.annotations['fragment_ends_at'] = int(fragend)
    fragment.annotations['SOURCE'] = YR_blast_qresult.target.split('/')[-1]
    fragment.annotations['contig'] = contig
    fragment.annotations['longID'] = contig + ':' + fragstart + ':' + fragend + ':' + YR_start
    
    #add YR features:
    YR_count = 1
    for hit in YR_blast_qresult:
        if hit.id == contig:
            for hsp in hit.hsps:
                if hsp_identifier == str(hsp.hit_start):
                    feature, YR_count = make_feature('YR', YR_blast_qresult, fragment ,hit, hsp, fragstart, YR_count)
                    fragment.features.append(feature)
    
    #add RT features:
    RT_count = 1
    for hit in RT_blast_qresult:
        if hit.id == contig + ':' + fragstart + ':' + fragend + ':' + YR_start:
            for hsp in hit.hsps:
                    feature, RT_count = make_feature('RT' , RT_blast_qresult, fragment ,hit, hsp, fragstart, RT_count)
                    fragment.features.append(feature)
   
    #add MT features:
    MT_count = 1
    for hit in MT_blast_qresult:
        if hit.id == contig + ':' + fragstart + ':' + fragend + ':' + YR_start:
            for hsp in hit.hsps:
                    feature, MT_count = make_feature('MT' , MT_blast_qresult, fragment ,hit, hsp, fragstart, MT_count)
                    fragment.features.append(feature)





#Searching repeats 





#Ugene#################################################################
def run_ugene(fragments_file, rep_direction,ident,maxdist, mindist ,minlength, outnamestart):
    d = 'false'
    if rep_direction == 'inverted':
        d = 'true'
    outfile = outnamestart + rep_direction + '.out'
    os.system('ugene find-repeats --in=' + fragments_file + '.fas --out=' + outfile + ' --identity=' + str(ident) + ' --max-distance=' + str(maxdist) + ' --inverted=' + d + ' --min-distance=' + str(mindist) + ' --min-length=' + str(minlength))
    return outfile
########################################################################




#print 'Searching direct repeats on fragments'
direct_repeats_file = run_ugene(fragments_file, 'direct',ident,maxdist, mindist ,minlength, outnamestart)
#print 'Searching inverted repeats on fragments'
inverted_repeats_file = run_ugene(fragments_file, 'inverted',ident,maxdist, mindist ,minlength, outnamestart)

#SeqIO does not parse ugene output. Therefore reading repeat features into SeqFeature objects.
#This will store all the seqfeature objects in format repeat_features[fragment][(seqfeature,seqfeature,...)]
#print 'Reading repeat features'
repeat_features = {}
repeat_name = 1
modifyers = {}
start_a = 0
end_a = 0
start_b = 0
end_b = 0
fragment = ""
direction = 'direct'
for filename in (direct_repeats_file, inverted_repeats_file):
    if filename == inverted_repeats_file:
        direction = 'inverted'
    lines = open(filename, 'r').readlines()
    readl = 0
    for line in lines:
        if line[0:9] == 'FASTA_HDR': #Get fragment name
            parts = re.split(r'\s+',line)
            fragment = ':'.join(parts[1].split(':')[0:-1])
            if not fragment in repeat_features.keys():
                repeat_features[fragment] = []
        elif line[0:8] == 'FEATURES': #Start to read feature lines
            readl = 1
        elif line[0:6] == 'ORIGIN': #Stop reading feature lines and ...
            readl = 0
            if not end_a == 0:      #If there were features, put the last one in repeat_features
                feature = SeqFeature(FeatureLocation(int(start_a)-1, int(end_a)), type="REPEAT", id = direction + '_' + str(repeat_name) + '.1') #Make a SeqFeature object for the first repeat partenr
                for mod in modifyers.keys():
                    feature.qualifiers[mod] = modifyers[mod] #add the qualifiers to the SeqFeature object
                feature.qualifiers['name'] = direction + '_' + str(repeat_name) + '.1'
                repeat_features[fragment].append(feature)
                feature = SeqFeature(FeatureLocation(int(start_b)-1, int(end_b)), type="REPEAT", id = direction + '_' + str(repeat_name) + '.2') #Make a SeqFeature object for the second repeat partenr
                for mod in modifyers.keys():
                    feature.qualifiers[mod] = modifyers[mod] #add the qualifiers to the SeqFeature object
                feature.qualifiers['name'] = direction + '_' + str(repeat_name) + '.2'
                #print ' Got repeat ' + feature.qualifiers['name'] + ' on fragment ' + fragment 
                repeat_features[fragment].append(feature)
                repeat_name = 1 #Roll back all parameters to null for the next fragment
                start_a = 0
                end_a = 0
                start_b = 0
                end_b = 0
                modifyers = {}
        #Now in between the strat and end for line reading:
        elif readl == 1:
            if line[5:16] == 'repeat_unit': 
                if not end_a == 0:    #Add previous repeat to a SeqFeature object in repeat_features
                    feature = SeqFeature(FeatureLocation(int(start_a)-1, int(end_a)), type="REPEAT", id = direction + '_' + str(repeat_name) + '.1')
                    for mod in modifyers.keys():
                        feature.qualifiers[mod] = modifyers[mod]
                    feature.qualifiers['name'] = direction + '_' + str(repeat_name) + '.1'
                    repeat_features[fragment].append(feature)
                    feature = SeqFeature(FeatureLocation(int(start_b)-1, int(end_b)), type="REPEAT", id = direction + '_' + str(repeat_name) + '.2')
                    for mod in modifyers.keys():
                        feature.qualifiers[mod] = modifyers[mod]
                    feature.qualifiers['name'] = direction + '_' + str(repeat_name) + '.2'
                    #print ' Got repeat ' + feature.qualifiers['name'] + ' on fragment ' + fragment
                    repeat_features[fragment].append(feature)
                    repeat_name += 1  
                    start_a = 0   #Roll back repeat values to null for next repeat
                    end_a = 0
                    start_b = 0
                    end_b = 0
                    modifyers = {}
                all_coords = line.split('(')[1]       #Read repeat coordinates
                all = re.sub('\)','',all_coords)
                all_coords = all.split(',')
                start_a = all_coords[0].split('..')[0]
                end_a = all_coords[0].split('..')[1]
                start_b = all_coords[1].split('..')[0]
                end_b = all_coords[1].split('..')[1]
                end_b = end_b.rstrip()
            elif line[21:22] == '/':                #Read repeat qualifiers
                line = line.rstrip('\n')
                modifyer = line.split('/')[1]
                key = modifyer.split('=')[0]
                value = modifyer.split('=')[1]
                value = re.sub('"','',value)
                modifyers[key] = value
#print repeat_features





#Insert repeats into SeqRecords
#print 'Annotating repeat features to fragments'
for fragment in fragment_stores:
    fragment_name = fragment.annotations['contig'] + ':' + str(fragment.annotations['fragment_starts_at']-1) + ':' + str(fragment.annotations['fragment_ends_at'])
    for key in repeat_features.keys():
        if fragment_name == key:
            for seqfeature in repeat_features[key]:
                fragment.features.append(seqfeature)
#Sort fragment.features based on feature.start
for fragment in fragment_stores:
    fragment.features.sort(key = lambda feature: int(feature.location.start))

                
#print 'Writing fragments to ' + fragments_file + '.embl'
#SeqIO.write(fragment_stores, fragments_file + '.embl', "embl")

    
#Make diagrams
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

#print 'Drawing'
for fragment in fragment_stores:
    diagram = GenomeDiagram.Diagram(fragment.description)
    track_for_features = diagram.new_track(1, name="Annotated Features")
    feature_set = track_for_features.new_set()
    color = ""
    sig = "BOX"
    lab = ""
    for feature in fragment.features:
        if feature.type == "DOMAIN":
            color = "gray"
            sig="BIGARROW"
            lab = feature.qualifiers['product']
        elif feature.type == "REPEAT":
            sig = "BOX"
            lab = feature.qualifiers['name']
            if feature.id[0:6] == 'direct':
                color = "green"
            else:
                color = "blue"
            
        feature_set.add_feature(feature, color=color , label=True, name=lab, sigil=sig, label_position="start",label_size=10, label_angle=90, arrowhead_length=10000)
    
    diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=7,start=0, end=len(fragment))
    diagram.write(work_path + 'fragments_diagrams/' + fragment.id + ".pdf", "PDF")





#Classify


import sys
import re
from cogent import LoadTable
from cogent.util.table import Table
column_headings = ['Fragment(species:contig:YR_start_on_contig)', 'Structure_Indication', 'start_on_contig', 'end_on_contig', 'Contig_length', 'Fragment_length', 'YR_evalue', 'RT_evalue', 'MT_evalue']
rows = []
TE_stores = []

for fragment in fragment_stores:
    YR_eval = '-'
    RT_eval = '-'
    MT_eval = '-'
    please_print = 0
    if fragment.id == '???':#test prints off
        please_print = 1   
    element = ''
    TE_related_features = []
    possibilities = []
    fragment_strand = 1
    YR_index = 0
    while YR_index < len(fragment.features):
        if fragment.features[YR_index].type == 'DOMAIN' and fragment.features[YR_index].qualifiers['product'] == 'YR':
            feature = fragment.features[YR_index]
            TE_related_features.append(feature)
            fragment_strand = feature.strand
            YR_eval = feature.qualifiers['evalue']
            break
        else:
            YR_index += 1
    if len(fragment.features) < 6: #minimal features in a YR retroelement (viper/Ngaro)
        possibilities = []
    else:
        
        if YR_index == 0 or YR_index == len(fragment.features) - 1:
            possibilities = []
        
        else:
        
            feature_b = fragment.features[YR_index + int(fragment_strand)]
            if please_print == 1:
                print feature_b
            feature_b_name = ''
            if feature_b.type == 'DOMAIN': #not expected because no domains are expected 3' of the YR
                possibilities = []
            elif feature_b.type == 'REPEAT':
                feature_b_name = feature_b.qualifiers['name']               #This is one of the TE repeat pairs
                                                                            #Check dirs1-like
                if 'inverted' in feature_b_name:
                    TE_related_features.append(feature_b)
                    possibilities = 'dirs1-like'
                    coordinates = []
                    coordinates.append(int(feature_b.location.start))
                    coordinates.append(int(feature_b.location.end))
                    feature_b_partners_name = feature_b_name[:-2] + '.2'
                    if '.2' in feature_b_name:
                        feature_b_partners_name = feature_b_name[:-2] + '.1'
                    for f in fragment.features:
                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_b_partners_name:
                            TE_related_features.append(f)
                            coordinates.append(int(f.location.start))
                            coordinates.append(int(f.location.end))
                            break
                    coordinates.sort()
                    if please_print == 1:
                        print coordinates
                    domains_nested_in_feature_b = 0
                    yr = 0
                    rt = 0
                    mt = 0
                    for f in fragment.features:
                        if f.type == 'DOMAIN' and coordinates[0] < int(f.location.start) and int(f.location.start) < coordinates[3]:
                            if f.qualifiers['product'] == 'YR':
                                yr = 1
                            elif f.qualifiers['product'] == 'RT':
                                rt = 1
                                TE_related_features.append(f)
                                RT_eval = f.qualifiers['evalue']
                            elif f.qualifiers['product'] == 'MT':
                                mt = 1
                                TE_related_features.append(f)
                                MT_eval = f.qualifiers['evalue']
                    domains_nested_in_feature_b = yr + rt + mt
                    if not domains_nested_in_feature_b == 3:                   #Checking that all the domains are nested in the first repeat pair as expected in dirs1
                        possibilities = []
                    feature_c_index = YR_index + 2 * (int(fragment_strand))
                    if feature_c_index < 0 or feature_c_index > len(fragment.features) - 1:
                        possibilities = []
                    else:
                        feature_c = fragment.features[feature_c_index]          #This is the second TE related TE pair
                        if please_print == 1:
                            print feature_c
                        if 'name' in feature_c.qualifiers.keys():
                            feature_c_name = feature_c.qualifiers['name']
                            coordinates = []
                            coordinates.append(int(feature_c.location.start))
                            coordinates.append(int(feature_c.location.end))
                            if not 'inverted' in feature_c_name:
                                possibilities = []
                            else:
                                feature_c_partners_name = feature_c_name[:-2] + '.2'
                                if '.2' in feature_c_name:
                                    feature_c_partners_name = feature_c_name[:-2] + '.1'
                                for f in fragment.features:
                                    if f.type == 'REPEAT' and f.qualifiers['name'] == feature_c_partners_name:
                                        coordinates.append(int(f.location.start))
                                        coordinates.append(int(f.location.end))
                                        break
                                coordinates.sort()
                                feature_c_good = 1
                                for f in TE_related_features:                              #making sure that none of the second repeat partners are nested in the first repeat
                                    if coordinates[1] < f.location.start < coordinates[2]: # allowing feature c and feature b (the repeats) to overlapp a little
                                        feature_c_good = 0
                                if feature_c_good == 0:
                                    possibilities = []
                                else:
                                    TE_related_features.append(feature_c)
                                    for f in fragment.features:
                                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_c_partners_name:
                                            TE_related_features.append(f)  
                                        
                elif 'direct' in feature_b_name:                                      #Now checking other DIRS with direct repeats
                    TE_related_features.append(feature_b)
                    possibilities = []
                    coordinates = []
                    coordinates.append(int(feature_b.location.start))
                    coordinates.append(int(feature_b.location.end))
                    feature_b_partners_name = feature_b_name[:-2] + '.2'
                    if '.2' in feature_b_name:
                        feature_b_partners_name = feature_b_name[:-2] + '.1'         
                    if please_print == 1:
                        print feature_b_partners_name
                    for f in fragment.features:
                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_b_partners_name:
                            TE_related_features.append(f)
                            coordinates.append(int(f.location.start))
                            coordinates.append(int(f.location.end))
                            break
                    coordinates.sort()
                    if please_print == 1:
                        print coordinates
                    yr = 0
                    rt = 0
                    mt = 0                                                             #Checking which domains are nested in first TE pair
                    for f in fragment.features:
                        if f.type == 'DOMAIN' and coordinates[0] < int(f.location.start) and int(f.location.start) < coordinates[3]:
                            if f.qualifiers['product'] == 'YR':
                                yr = 1
                            elif f.qualifiers['product'] == 'RT':
                                rt = 1
                                TE_related_features.append(f)
                                RT_eval = f.qualifiers['evalue']
                            elif f.qualifiers['product'] == 'MT':
                                mt = 1
                                TE_related_features.append(f)
                                MT_eval = f.qualifiers['evalue']
                    if yr == 1 and rt == 0 and mt == 0:                              # only YR nested in first repeat, pat expected
                        feature_c_index = YR_index + (-1) * (int(fragment_strand))
                        if 0 <= feature_c_index and feature_c_index < len(fragment.features) - 1:
                            feature_c = fragment.features[feature_c_index]
                            TE_related_features.append(feature_c)
                            if please_print == 1:
                                print feature_c
                            if feature_c.type == 'REPEAT':
                                feature_c_name = feature_c.qualifiers['name']
                                coordinates = []
                                coordinates.append(int(feature_c.location.start))
                                coordinates.append(int(feature_c.location.end))
                                if not 'direct' in feature_c_name:
                                    possibilities = []
                                else:
                                    feature_c_partners_name = feature_c_name[:-2] + '.2'
                                    if '.2' in feature_c_name:
                                        feature_c_partners_name = feature_c_name[:-2] + '.1'
                                    for f in fragment.features:
                                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_c_partners_name:
                                            coordinates.append(int(f.location.start))
                                            coordinates.append(int(f.location.end))
                                            TE_related_features.append(f)
                                            break
                                    coordinates.sort()
                                    if please_print == 1:
                                        print coordinates
                                    nested_repeat_partners = 0
                                    if please_print == 1:
                                        print TE_related_features
                                    for fb in TE_related_features:
                                        if fb.type == 'REPEAT':
                                            if coordinates[0] < fb.location.start and fb.location.end  < coordinates[3]:   #Making sure that exectly one partner of first repeat pair is nested in the second
                                                nested_repeat_partners += 1
                                    if nested_repeat_partners == 1:
                                        if please_print == 1:
                                            print 'ok'
                                        for f in TE_related_features:
                                            if f.type == 'DOMAIN' and f.qualifiers['product'] == 'YR':                     #making sure that YR not in second repeat and MT and RT are.       
                                                if please_print == 1:
                                                    print 'ok.1'
                                                if not coordinates[0] < f.location.start < coordinates[3]:                 
                                                    if please_print == 1:
                                                        print 'first ok'
                                                    for frt in fragment.features:
                                                        if frt.type == 'DOMAIN' and frt.qualifiers['product'] == 'RT':
                                                            if coordinates[0] < frt.location.start and frt.location.start < coordinates[3]:
                                                                TE_related_features.append(frt)
                                                                RT_eval = frt.qualifiers['evalue']
                                                                if please_print == 1:
                                                                    print 'second ok'
                                                                for fmt in fragment.features:
                                                                    if fmt.type == 'DOMAIN' and fmt.qualifiers['product'] == 'MT':
                                                                        if please_print == 1:
                                                                            print fmt
                                                                        if coordinates[0] < fmt.location.start and fmt.location.start < coordinates[3]:
                                                                            TE_related_features.append(fmt)
                                                                            MT_eval = fmt.qualifiers['evalue']
                                                                            possibilities = 'pat1-like'
                                                

                            
                    elif yr == 0 and rt == 0 and mt == 0:                         #No domains nested in first repeat, toc expected
                        feature_c_index = YR_index + 2 * (int(fragment_strand))
                        if 0 <= feature_c_index < len(fragment.features) - 1:
                            feature_c = fragment.features[feature_c_index]
                            if feature_c.type == 'REPEAT':
                                feature_c_name = feature_c.qualifiers['name']
                                coordinates = []
                                coordinates.append(int(feature_c.location.start))
                                coordinates.append(int(feature_c.location.end))
                                TE_related_features.append(feature_c)
                                if not 'direct' in feature_c_name:
                                    possibilities = []
                                else:
                                    feature_c_partners_name = feature_c_name[:-2] + '.2'        #Cecking that only one partner of first repeat is in second repeat
                                    if '.2' in feature_c_name:
                                        feature_c_partners_name = feature_c_name[:-2] + '.1'
                                    for f in fragment.features:
                                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_c_partners_name:
                                            coordinates.append(int(f.location.start))
                                            coordinates.append(int(f.location.end))
                                            TE_related_features.append(f)
                                            break
                                    coordinates.sort()
                                    nested_repeat_partners = 0
                                    for fb in TE_related_features:
                                        if fb.type == 'REPEAT':
                                            if coordinates[0] < fb.location.start and fb.location.end < coordinates[3]:
                                                nested_repeat_partners += 1
                                    if nested_repeat_partners == 1:                               #making sure that all domains are in second repeat pair
                                        for f in TE_related_features:
                                            if f.type == 'DOMAIN' and f.qualifiers['product'] == 'YR':
                                                if (coordinates[0] < f.location.start < coordinates[3]):
                                                    for frt in fragment.features:
                                                        if frt.type == 'DOMAIN' and frt.qualifiers['product'] == 'RT':
                                                            if coordinates[0] < frt.location.start < coordinates[3]:
                                                                TE_related_features.append(frt)
                                                                RT_eval = frt.qualifiers['evalue']
                                                                found_fmt = 0
                                                                for fmt in fragment.features: 
                                                                    if fmt.type == 'DOMAIN' and fmt.qualifiers['product'] == 'MT':
                                                                        found_fmt = 1
                                                                        if coordinates[0] < fmt.location.start < coordinates[3]:
                                                                            TE_related_features.append(fmt)
                                                                            MT_eval = fmt.qualifiers['evalue']
                                                                            possibilities = 'PAT'
                                                                    if found_fmt == 0:
                                                                        possibilities = 'viper_or_Ngaro'
                    elif yr == 0 and rt == 1 and mt == 1:                            #Cheking kangaroo
                        feature_c_index = YR_index + 2* (int(fragment_strand))
                        if 0 <= feature_c_index < len(fragment.features) - 1:
                            feature_c = fragment.features[feature_c_index]
                            if feature_c.type == 'REPEAT':
                                feature_c_name = feature_c.qualifiers['name']
                                coordinates = []
                                coordinates.append(int(feature_c.location.start))
                                coordinates.append(int(feature_c.location.end))
                                TE_related_features.append(feature_c)
                                if not 'direct' in feature_c_name:
                                    possibilities = []
                                else:
                                    feature_c_partners_name = feature_c_name[:-2] + '.2'
                                    if '.2' in feature_c_name:
                                        feature_c_partners_name = feature_c_name[:-2] + '.1'
                                    for f in fragment.features:
                                        if f.type == 'REPEAT' and f.qualifiers['name'] == feature_c_partners_name:
                                            coordinates.append(int(f.location.start))
                                            coordinates.append(int(f.location.end))
                                            TE_related_features.append(f)
                                            break
                                    coordinates.sort()
                                    nested_repeat_partners = 0
                                    for fb in TE_related_features:
                                        if fb.type == 'REPEAT':
                                            if coordinates[0] < fb.location.start and fb.location.end < coordinates[3]:
                                                nested_repeat_partners += 1
                                    if nested_repeat_partners == 1:
                                        for f in TE_related_features:
                                            if f.type == 'DOMAIN' and f.qualifiers['product'] == 'YR':
                                                if (coordinates[0] < f.location.start < coordinates[3]):
                                                    for frt in fragment.features:
                                                        if frt.type == 'DOMAIN' and frt.qualifiers['product'] == 'RT':
                                                            if not coordinates[0] < frt.location.start < coordinates[3]:
                                                                for fmt in fragment.features:
                                                                    if fmt.type == 'DOMAIN' and fmt.qualifiers['product'] == 'MT':
                                                                        if not coordinates[0] < fmt.location.start < coordinates[3]:
                                                                            possibilities = 'kangaroo'
 

    if len(possibilities) == 0:               #If non of the above, collecting some info, eg, number of domains and their strand
                        yr = 0
                        yr_strand = 1
                        rt = 0
                        rt_strand = 1
                        mt = 0
                        mt_strand = 1
                        for feature in fragment.features:
                            if feature.type == 'DOMAIN':
                                if feature.qualifiers['product'] == 'YR':
                                    yr = 1
                                    yr_strand = feature.strand
                                elif feature.qualifiers['product'] == 'RT':
                                    rt = 1
                                    rt_strand = feature.strand
                                    TE_related_features.append(feature)
                                    RT_eval = feature.qualifiers['evalue']
                                elif feature.qualifiers['product'] == 'MT':
                                    mt = 1
                                    mt_strand = feature.strand
                                    TE_related_features.append(feature)
                                    MT_eval = feature.qualifiers['evalue']
                        if yr + rt + mt == 3:
                            if yr_strand == rt_strand == mt_strand:
                                possibilities = ('DIRS')
                            else:
                                possibilities = ('domains_inverted')
                        elif yr + rt + mt == 2:
                            possibilities = ('two_domains')
                        elif yr + rt + mt == 1:
                            possibilities = ('YR_only')
                       
    #Draw element diagram
    diagram = GenomeDiagram.Diagram(fragment.description + '_TE')
    track_for_features = diagram.new_track(1, name="Annotated Features")
    feature_set = track_for_features.new_set()
    color = ""
    sig = "BOX"
    lab = ""
    for feature in TE_related_features:
            if feature.type == "DOMAIN":
                color = "gray"
                sig="BIGARROW"
                lab = feature.qualifiers['product']
            elif feature.type == "REPEAT":
                sig = "BOX"
                lab = feature.qualifiers['name']
                if feature.id[0:6] == 'direct':
                    color = "green"
                else:
                    color = "blue"
            
            feature_set.add_feature(feature, color=color , label=True, name=lab, sigil=sig, label_position="start",label_size=10, label_angle=90, arrowhead_length=10000)
    
    diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=7,start=0, end=len(fragment))
    el = re.sub(' ','_',element)
    diagram.write(work_path + 'TEs_diagrams/' +  fragment.id + "_TE.pdf", "PDF")                             
                    
    #determine element start and end of element and put in fragment SeqRecord object                    
    all_te_coords = []
    
    for f in TE_related_features:        # 23.9.13 Now I only include domains that are in between te related repeats. This was done to make sure that I dont take RTs that dont belong to the TRE (= thyrosine recombinase element)
        if f.type == "REPEAT":
            all_te_coords.append(f.location.start)
            all_te_coords.append(f.location.end)
        elif f.type == "DOMAIN" and f.qualifiers['product'] == 'YR':
            all_te_coords.append(f.location.start)
            all_te_coords.append(f.location.end)
    all_te_coords.sort()
    start = all_te_coords[0]
    end = all_te_coords.pop()
    TE = SeqFeature(FeatureLocation(start, end), type="TE")
    TE.qualifiers['element'] = possibilities
    fragment.features.append(TE)
    fragment_length = fragment.annotations['fragment_ends_at'] - fragment.annotations['fragment_starts_at'] + 1
            
    rows.append([genome_code + ':' + str(fragment.annotations['contig']) + ':' + fragment.id.split('_')[1], possibilities, fragment.annotations['fragment_starts_at'], fragment.annotations['fragment_ends_at'], len(seq_for_contig(fragment.annotations['contig'] ,recs)), fragment_length, YR_eval, RT_eval, MT_eval])
    fragment.annotations['element'] = possibilities
    

        
                    
                    
                    
# Print results table
def formatcol(value):
    if isinstance(value, float):
        val = "%.4e" % value
    else:
        val = str(value)
    return val

t = Table(header = column_headings, rows = rows, column_templates = dict(YR_evalue=formatcol, RT_evalue=formatcol, MT_evalue=formatcol))
t.Title = "YR element hits in " + genome_code
t.Legend = 'Indications are structural and not phylogenetic.\nIndications depend on Ugene parameters.\n1. pat1-like: ==>=RT=MT=>>=>=YR=>>=.\n2. kangaroo: ==>=RT=MT=>>=>=RY=>>=.\n3. PAT: =>RT=MT=YR=>>=>=>>.\n4. dirs1-like: =>=RT=MT=YR=<>>=<<=.\n5. Ngaro or Viper: =>RT=YR=>>=>=>>=.'
print t
t.writeToFile(outnamestart + 'table.out')
                    



#Write protein fasta

#If multiple RT or multiple MT on a single fragment, if the distance between them is less than the allowed intron length and both of them are within the TE coordinates
#and they are on the same strand, they will be concatenated in the fasta file
#If they are both in the element but the distance between them is larger than allowed intron they will be written to the fasta file as individual sequences and only the 
#sequence marked as 'serial_on_frag=1' will be concatenated with the other domains
#a warning will be made and it is important to check manually because another solution may be better, eg, choose number 2 or concatenate anyways.


for I in ('YR','RT','MT'):
    protfasta = open(outnamestart + I + '_aa.fas', 'wt')
    for fragment in fragment_stores:
        te_start = ''
        te_end = ''
        for feature in fragment.features:
            if feature.type == 'TE':
                te_start = feature.location.start
                te_end = feature.location.end
        Is = []
        for feature in fragment.features:
            if 'product' in feature.qualifiers.keys() and te_start <= feature.location.start < te_end:
                if I in feature.qualifiers['product']:
                    Is.append(feature)
        if len(Is) > 1:
            sequence = ''
            coordinates = []
            for i in Is:
                coordinates.append(i.location.start)
                coordinates.append(i.location.end)
            coordinates.sort()
            if coordinates[2] - coordinates[1] < intron:
                if Is[0].strand == Is[1].strand == 1:
                    sequence = Is[0].qualifiers['translation'] + Is[1].qualifiers['translation']
                elif Is[0].strand == Is[1].strand == -1:
                    sequence = Is[1].qualifiers['translation'] + Is[0].qualifiers['translation']
            if len(sequence) > 1:
                protfasta.write('>' + genome_code + ':' + fragment.annotations['contig'] + ':' + fragment.id.split('_')[1] + ':' + '1' + ':' + fragment.annotations['element'] + '\n' + str(sequence) + '\n')
                print 'WARNING: Two ' + I + ' hits were concatenated in ' + genome_code + ':' + fragment.annotations['contig'] + ':' + fragment.id.split('_')[1]
            else:
                for i in Is:
                    protfasta.write('>' + genome_code + ':' + fragment.annotations['contig'] + ':' + fragment.id.split('_')[1] + ':' + str(i.qualifiers['serial_on_frag']) + ':' + fragment.annotations['element'] + '\n' + str(i.qualifiers['translation']) + '\n')
                    if i.qualifiers['serial_on_frag'] > 1:
                        print 'WARNING: Multiple ' + I + ' hits on ' + genome_code + ':' + fragment.annotations['contig'] + ':' + fragment.id.split('_')[1] + ' are over ' + str(intron) + ' bp apart and will not be concatenated. Only serial_on_frag=1 will be used in tree.'
        elif len(Is) > 0:
            protfasta.write('>' + genome_code + ':' + fragment.annotations['contig'] + ':' + fragment.id.split('_')[1] + ':' + str(Is[0].qualifiers['serial_on_frag']) + ':' + fragment.annotations['element'] + '\n' + str(Is[0].qualifiers['translation']) + '\n')
    protfasta.close()

SeqIO.write(fragment_stores, fragments_file + '.embl', "embl")    
   
