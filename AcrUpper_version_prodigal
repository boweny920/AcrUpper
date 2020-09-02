###Useage: python3 Pipeline_1.py GCF_Folder OutputFolder[> log_file]

## Must have Known_Acr.dmnd file and "hmmpress hmmdb_HTH_domain" output file and the assembly_summary file present in the directory of this script!!
## This script will find Acr homologs with diamond, Find potential Acr-Aca loci and find perform a hmmscan with the 65 published HTH domain againest those loci.
## Output of the script: Acr homologs fasta file; Potential Acr-Aca loci file; Candidate Aca fasta file; Diamond output file hmmscan output file
## look for file "!!!!! Potential Acr-Aca loci found and saved in %s"; "$$$$$ Aca candidate found and saved in file %s" for final output results


import subprocess
import sys
from Bio import SeqIO
import os
import pprint as pp

##Filteration Functions:
def make_output_folder(GCF_Name):
    os.makedirs(GCF_Name)
    return os.path.dirname(GCF_Name)

def find_faa(GCF_directory):
    ##modified for prodigal
    faa = subprocess.run(["find", GCF_directory, "-maxdepth", "15", "-name", "*fna_prodigal_out.faa"],stdout=subprocess.PIPE).stdout.decode('utf-8')
    faa = str(faa).strip()
    return faa

def find_gff(GCF_directory):
    ##modified for prodigal
    gff = subprocess.run(["find", GCF_directory, "-maxdepth", "15", "-name", "*fna_prodigal_out.gff"],stdout=subprocess.PIPE).stdout.decode('utf-8')
    gff = str(gff).strip()
    return gff

def get_file_directory_path(file_path):
    directory_path=str(file_path).rsplit("/",1)
    return directory_path[0]

def get_faa_name(file_path):
    directory_path = str(file_path).rsplit("/", 1)
    return directory_path[1]

def length(list):
    gene_length=int(list[4])-int(list[3]) + 1
    return gene_length

def contig(list):
    contig=list[0]
    return contig

def strand(list):
    gene_strand=list[6]
    return gene_strand

def proteinInfo(list):
    ##modified to fit prodigal output
    return [list[0], list[-1]]

def proteinID(IDinfo_prodiagal):
    ##modified to fit prodigal output
    IDinfo = IDinfo_prodiagal[-1].split(";")
    for v in IDinfo:
        if "ID=" in v:
            ID = str(v).split("_",1)
            ProID=str(IDinfo_prodiagal[0])+"_"+ID[-1]
            return ProID

def proteinInfo_List_process(proteinInfo_list):
    combined_list=[]
    for sublist in proteinInfo_list:
        loci_list_proteins_WithPseudo = []
        for number,value in enumerate(sublist):
            if "pseudo=true" not in value:
                loci_list_proteins_WithPseudo.append([int(number),proteinID(value)])
            elif "pseudo=true" in value:
                loci_list_proteins_WithPseudo.append([int(number),"pseudo"])
        pseudo_gene_number=0
        for x in loci_list_proteins_WithPseudo:pseudo_gene_number=pseudo_gene_number+x.count("pseudo")
        if len(loci_list_proteins_WithPseudo) - pseudo_gene_number >= 2 :
        ##see how many genes left without the pseudo gene
            combined_list.append(loci_list_proteins_WithPseudo)
    return combined_list

def find_GCF(GCF_directory):
    faa_name=get_faa_name(find_faa(GCF_directory))
    name=faa_name.split("_",2)
    GCF=name[0]+"_"+name[1]
    return GCF

def species(assembly_dic,GCF):
    if GCF in assembly_dic.keys():
        return assembly_dic[GCF]
    if GCF not in assembly_dic.keys():
        print("%s not in assembly_summary"%GCF)

def result_check_list(list,GCF,assembly_dic):
    gene_strand=strand(list)
    gene_length=length(list)
    species_name=species(assembly_dic,GCF)
    if "pseudo=true" not in list[-1]:ProID=proteinID(proteinInfo(list))
    elif "pseudo=true" in list[-1]:ProID="pseudo"
    # output result in format ProteinID; contig;gene strand; gene length；gene start;gene end, in a list
    check_result = [ProID, contig(list), gene_strand, gene_length, list[3], list[4],species_name]
    return check_result

def loci_select_before_diamond(file_dic,GCF,assembly_dic):
    ###################################################################
    ###   Get the regions that contain at least two genes, each gene must
    ###   be on the same strand, each gene must be less than 200 aa or 600 nucleotides
    key=0
    lst=[]
    loci_list=[]
    lst_result_check=[]
    loci_list_result_check=[]
    while key in file_dic.keys():
        ##First filter gene length:
        if length(file_dic[key]) < 600:
            #This is working!!
            if key-1 in file_dic.keys() and strand(file_dic[key]) == strand(file_dic[key-1]) and contig(file_dic[key]) == contig(file_dic[key-1]) and length(file_dic[key-1]) < 600:
                #print(contig(file_dic[key]))
                lst.append(proteinInfo(file_dic[key]))
                # This is for result checking
                lst_result_check.append(result_check_list(file_dic[key],GCF,assembly_dic))
                key=key+1
            elif key+1 in file_dic.keys() and strand(file_dic[key])==strand(file_dic[key+1]) and contig(file_dic[key]) == contig(file_dic[key+1]) and length(file_dic[key+1]) < 600:
                if len(lst) >= 2:
                    #if previous lst has 2 or more genes, append the lst to the loci_list
                    loci_list.append(lst)
                    loci_list_result_check.append(lst_result_check)
                lst = []
                lst_result_check = []
                lst.append(proteinInfo(file_dic[key]))
                lst.append(proteinInfo(file_dic[key+1]))
                #This is for result checking
                lst_result_check.append(result_check_list(file_dic[key],GCF,assembly_dic))
                lst_result_check.append(result_check_list(file_dic[key+1],GCF,assembly_dic))

                key=key+2
            else:
                key = key + 1
                if len(lst) >= 2:
                    # print(lst)
                    loci_list.append(lst)
                    loci_list_result_check.append(lst_result_check)
                lst = []
                lst_result_check = []
        else:
            key=key+1
            if len(lst) >= 2:
                # print(lst)
                loci_list.append(lst)
                # This is for result checking
                loci_list_result_check.append(lst_result_check)
            lst=[]
            lst_result_check = []
    else:
        if len(lst) >= 2:
            # print(lst)
            loci_list.append(lst)
            # This is for result checking
            loci_list_result_check.append(lst_result_check)
        lst = []
        lst_result_check = []

    loci_list_w_pseudo=proteinInfo_List_process(loci_list)
    # pp.pprint(loci_list_w_pseudo)
    # pp.pprint(loci_list_result_check)
    return loci_list_w_pseudo,loci_list_result_check

def is_non_zero_file(fpath):
    #Check if file is empty or does not exsit
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def make_and_check_output_directory(dpath):
    if os.path.isdir(dpath) is False:
        os.makedirs(dpath)
    return str(dpath)

def run_diamond(query_file,outputfile):
    ##Use this as dmond_out
    subprocess.Popen(['diamond','blastp','-d','Known_Acr','-f','6','--quiet','--more-sensitive','-q',query_file,'-p','5','-o',str(outputfile)+'.dmnd_out']).wait()
    ## Output format 6 columns: Default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    return str(outputfile)+'.dmnd_out'

def parse_diamond_get_proteinID(dmond_out):
    ##Use with run_diamond; Use on NON-EMPTY diamond output file
    ##Return also a dict of Protein ID to Acr type
    dic_acr = {}
    protein_ID_list = []
    Pro_ID = subprocess.Popen("awk '{print $1,$2}' %s" % dmond_out, shell=True, stdout=subprocess.PIPE)
    for line in Pro_ID.stdout:
        line = line.rstrip().decode('utf-8').split()
        proID = str(line[0])
        protein_ID_list.append(proID)
        dic_acr.setdefault(proID, str(line[1]))
    return protein_ID_list, dic_acr

def faa_file_wrote(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    ##The one_acr_aca_locus is a list of proteinID.
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory,"w") as newfile:
        for proID in one_acr_aca_locus:
            SeqIO.write(record_dict[proID],newfile,"fasta")

def faa_file_wrote_potential_AcrAca(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory, "w") as newfile:
        for gene_order,proID in one_acr_aca_locus:
            if proID in record_dict.keys():
                SeqIO.write(record_dict[proID], newfile, "fasta")

def faa_file_wrote_diamond_special(one_acr_aca_locus,faa_file,newfile_name_dirctory,Acr_protein_NP_list):
    ##This function takes the proteins that are not Acrs and write them into a new faa file
    ##The one_acr_aca_locus is a list of [gene order, proteinID].
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with gene orders with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    unhmmer_aca = []
    with open(newfile_name_dirctory, "w") as newfile:
        acr_location=[]
        for order, proID in one_acr_aca_locus:
            if proID in Acr_protein_NP_list: acr_location.append(int(order)+1)
            # "+1" is to get rid of 0, which makes the stops the function when occured
        for order1, proID1 in one_acr_aca_locus:
            order1=order1+1
            # "+1" is to get rid of 0, which makes the stops the function when occured
            #if proID1 not in Acr_protein_NP_list and any(v for v in acr_location if abs(v-int(order1))<=4):
            if any(v for v in acr_location if abs(v - int(order1)) <= 4):
                if proID1 in record_dict.keys():
                    SeqIO.write(record_dict[proID1], newfile, "fasta")
                    unhmmer_aca.append(proID1)
                ##Check the gene distance between Acr and Aca, it should be less than 3 genes
    return unhmmer_aca

def run_hmmscan(faafile):
    hmm_outfile=str(faafile)+".hmmout"
    subprocess.Popen(['hmmscan','--domtblout',hmm_outfile,'-o','log.hmm','--noali','-E','1e-2','hmmdb_HTH_domain',faafile]).wait()
    return hmm_outfile

def parse_hmmOutfile(hmm_outfile):
    parsed_hmm_outfile = str(hmm_outfile) + ".Coverage_parsed"
    subprocess.Popen(["grep -v '#' %s | awk '($19-$18)/($17-$16)>0.5 {print $0}' > %s" % (hmm_outfile, parsed_hmm_outfile)],shell=True).wait()
    return parsed_hmm_outfile

def potential_new_aca_faa_filemake(parsed_hmm_outfile,newAcrAca_faaFile,unhmmer_aca):
    ## Make a fasta file of potential Aca amino acid sequences
    ## Make a list of Aca, for checking results
    Aca_Pro_lst=[]
    record_dict = SeqIO.to_dict(SeqIO.parse(newAcrAca_faaFile, "fasta"))
    newfile_name= parsed_hmm_outfile+".new_found_ACA.faa"
    with open(newfile_name,"w") as new:
        acaout=subprocess.Popen(["awk '{print $4}' %s |sort -u "%parsed_hmm_outfile],shell=True, stdout=subprocess.PIPE)
        for ID in acaout.stdout:
            ID=ID.strip().decode('utf-8')
            if ID in unhmmer_aca:
                Aca_Pro_lst.append(ID)
                # If Aca is not in unhmmer, then it should not be considered as a potential Aca(here is to make sure that condiction "gene distance < 3" is applied)
                SeqIO.write(record_dict[ID], new, "fasta")
                #This is to consider the < than 3 genes condition, which had been implemented before in "faa_file_wrote_diamond_special" function
    return newfile_name,Aca_Pro_lst

def final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,newfile_name,GCF):
    #generation of final result checking file
    #Table will be tsv format
    #### GCF_number    ProteinID   Contig  Strand  Length  Start   End     Species     Acr     Aca ####
    Acr_Aca_loci_by_GBA = []
    for i in Acr_homolog_candidate_list_resultCheck:
        if any(v for v in i if v[0] in Aca_Pro_lst):
            for v in i:
                # Check which ProteinID is the Aca homolog
                if v[0] in Aca_Pro_lst:
                    v.append("Aca_protein")
                else:
                    v.append("NA")
            Acr_Aca_loci_by_GBA.append(i)
    # pp.pprint(Acr_Aca_loci_by_GBA)
    # Acr_Aca_loci_checkResult = output_1_DirPath + "/AcrAca_loci_ORDER=" + str(order) + ".check_Result"
    with open(newfile_name, "w") as newfile:
        for pro in Acr_Aca_loci_by_GBA:
            for value in pro:
                newfile.write("%s\t"%str(GCF))
                for info in value: newfile.write("%s\t" % str(info))
                newfile.write("\n")
##########################################################################################

GCF_directory=sys.argv[1]
##This is the GCF file directory
output_folder=make_and_check_output_directory(str(sys.argv[2]))
##Specifying the out put folder

faa=find_faa(GCF_directory)
gff=find_gff(GCF_directory)

### Write the assembly information into a dictionary
assembly_dic={}
assem=subprocess.Popen(["awk 'NR>2 {print $0}' assembly_summary.txt.Refseq_bacteriaGenomes_10-21-19"],shell=True,stdout=subprocess.PIPE)
for line in assem.stdout:
    line=line.strip().decode('utf-8')
    line=line.split("\t")
    assembly_dic.setdefault(line[0],line[7])

GCF=find_GCF(GCF_directory)

if faa is not None and gff is not None:
##This is the protein faa file
    GCF_path=get_file_directory_path(faa)
    faa_name=get_faa_name(faa)
    sub_outputfolder_path=output_folder+"/"+str(faa_name)
    os.makedirs(sub_outputfolder_path)

    file_list=[]
    cmd=['grep "CDS" %s'%gff]
    file=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    for line in file.stdout:
        line = line.decode('utf-8').rstrip().split('\t')
        file_list.append(line)
    file.wait()
    file_dic={}
    for i in enumerate(file_list):
        file_dic.setdefault(i[0],i[1])

##The GCF and assembly_dic is for result checking, specifically for making the species names
    loci_list,loci_list_result_check=loci_select_before_diamond(file_dic,GCF,assembly_dic)
    diamond_outfile=run_diamond(faa,sub_outputfolder_path+"/"+str(faa_name))

    if is_non_zero_file(diamond_outfile) is not False:
        Acr_homolog_candidate_list=[]
        #For result checking
        Acr_homolog_candidate_list_resultCheck = []
        protein_NP_list, dic_Acr = parse_diamond_get_proteinID(diamond_outfile)
        ##Write the Acr homologs to a fasta file
        faa_file_wrote(protein_NP_list,faa,sub_outputfolder_path+"/"+str(faa_name)+".Acr_homologs.faa")
        for i in loci_list:
            # Note: v[1] would be the protein ID not the gene order
            if any(v for v in i if v[1] in protein_NP_list):Acr_homolog_candidate_list.append(i)
            ##To see if any sub-list of loci_list has values in protein_NP_list.

        ### For result checking, v[0] is the Protein ID
        for i in loci_list_result_check:
            if any(v for v in i if v[0] in protein_NP_list):
                for v in i:
                    #Check which ProteinID is the Acr homolog
                    if v[0] in protein_NP_list:v.append(dic_Acr[v[0]])
                    else:v.append("NA")
                Acr_homolog_candidate_list_resultCheck.append(i)
        ####
        if len(Acr_homolog_candidate_list) > 0:
            #Put the new faa files in a directory called "Acr_homolog_Acr-Aca_candidate"
            output_1_DirPath = sub_outputfolder_path + "/Acr_homolog_Acr-Aca_candidate"
            os.makedirs(output_1_DirPath)
            #Write each locus with Acr homolog into individual fasta files
            for order, locus in enumerate(Acr_homolog_candidate_list):
                potential_Acr_Aca_filename=output_1_DirPath+"/Acr_Aca_candidate_with_Acr-Homolog_ORDER="+str(order)+".faa"
                potential_Aca_filename = output_1_DirPath + "/Aca_candidate_in_Acr_Aca_candidate_with_Acr-Homolog_ORDER=" + str(order) + ".faa"
                unhmmer_aca=faa_file_wrote_diamond_special(locus,faa,potential_Aca_filename,protein_NP_list)
                faa_file_wrote_potential_AcrAca(locus, faa, potential_Acr_Aca_filename)
                print("!!!!! Potential Acr-Aca loci found and saved in %s"%potential_Acr_Aca_filename)
                if is_non_zero_file(potential_Aca_filename):
                    print("!!!!! Potential Aca genes found and saved in %s" % potential_Aca_filename)
                ##run hmmscan on this "potential_Aca_Aca_filename"
                hmm_outfile=run_hmmscan(potential_Acr_Aca_filename)
                if is_non_zero_file(hmm_outfile) is not False:
                    ##parse the hmm_outfile based on coverage>0.5
                    parsed_hmm_outfile=parse_hmmOutfile(hmm_outfile)
                    if is_non_zero_file(parsed_hmm_outfile) is not False:
                        #Return the potential Aca sequences in a fasta file
                        # aca_candidate_file=potential_new_aca_faa_filemake(parsed_hmm_outfile, potential_Acr_Aca_filename)
                        ##This code generate Traceback (most recent call last):
                            #   File "Pipeline_1_V3.py", line 258, in <module>
                            #     aca_candidate_file=potential_new_aca_faa_filemake(parsed_hmm_outfile, potential_Acr_Aca_filename)
                            #   File "Pipeline_1_V3.py", line 186, in potential_new_aca_faa_filemake
                            #     record_dict = SeqIO.to_dict(SeqIO.parse(newAcrAca_faaFile, "fasta"))
                            #   File "/home/yangbowen/.local/lib/python3.6/site-packages/Bio/SeqIO/__init__.py", line 813, in to_dict
                            #     raise ValueError("Duplicate key '%s'" % key)
                            # ValueError: Duplicate key 'WP_000956747.1'
                        #The below code may solve this problem
                        #The below code also generated a Aca_Pro list for checking results
                        aca_candidate_file,Aca_Pro_lst = potential_new_aca_faa_filemake(parsed_hmm_outfile,faa,unhmmer_aca)
                        if is_non_zero_file(aca_candidate_file) is not False:
                            print("$$$$$ Aca candidate found and saved in file %s"%aca_candidate_file)

                        ####For result checking#####
                        Acr_Aca_loci_checkResult = output_1_DirPath + "/AcrAca_loci_ORDER=" + str(order) + ".check_Result"
                        final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,Acr_Aca_loci_checkResult,GCF)
                        #####For result checking#####
        elif len(Acr_homolog_candidate_list) is 0:
            print("##### %s have %s number of Acr homologs, but no qualified Acr-Aca loci #####"%(faa_name,len(protein_NP_list)))
    else:print("##### No Acr homologs in annotated genome '%s' #####"%faa_name)

else: print("##### %s has no gff or no faa file #####"%GCF_directory)
