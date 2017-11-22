import yaml

yaml_file = "_global_projectb_scratch_andreopo_ACDC_acdc_build_SAG_contam_Meiothermus_ruber_SAG.Pedobacter_heparinus_SAG_Meiothermus.Pedobacter.fastq.CAREFUL_scaffolds.fasta.min2kb.fasta.yaml"

with open (yaml_file, 'r') as f:
    doc=yaml.load(f)
    # use safe_load instead load
    # dataMap = yaml.safe_load(f)


method=doc['cluster_estimates']['most_likely_clustering']['method']

num_clusters = doc['cluster_estimates']['most_likely_clustering']['estimated_k']

contig_names= doc['fasta_stats']['included_contigs']
contigs_to_unsupclusters = {}
unsupclusters_to_contigs = {}

unsup_clusters = []
all_k_clusterings = doc['cluster_estimates'][method]['assignments']###['assignment'][num_clusters]
for a in all_k_clusterings:
   if a.get('assignment').get('k') == num_clusters:
       unsup_clusters = a.get('assignment').get('labels')
count = 0
for cluster in unsup_clusters:
   contig_name = contig_names[count]
   contigs_to_unsupclusters[contig_name] = cluster
   if not cluster in unsupclusters_to_contigs:
      unsupclusters_to_contigs[cluster] = []
   unsupclusters_to_contigs[cluster].append(contig_name)
   count+=1

cluster_to_max_confidence = {}
cluster_to_class = {}
taxa_contigs =  doc['taxonomies']  ###a dictionary with the contig names as keys
for t in taxa_contigs:
   cluster = contigs_to_unsupclusters.get(t)
   conf = taxa_contigs.get(t).get("confidence")
   classif = taxa_contigs.get(t).get("state")
   if not cluster in cluster_to_class or conf > cluster_to_max_confidence[cluster]:
     cluster_to_class[cluster] = classif
     cluster_to_max_confidence[cluster] = conf

for c in cluster_to_class:
   classif = cluster_to_class.get(c)
   file_classif = open(classif + ".names", 'a')
   contigs = unsupclusters_to_contigs.get(c)
   for name in contigs:
      file_classif.write(name + '\n')
   file_classif.close()


###for i in `ls *.names` ; do filterbyname.sh names=$i in=$fasta out=$i.fasta ; done
