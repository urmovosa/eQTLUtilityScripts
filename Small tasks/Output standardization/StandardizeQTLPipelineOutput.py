import cProfile
import gzip
import argparse

# Set up the usage of external arguments:
parser = argparse.ArgumentParser(description = 'Update eQTL output files to new cohort names and order.')

parser.add_argument('--mapping_file', metavar = 'name mapping file',
                    help = 'File with the new cohort name. Three columns: name_old- cohort name in the original eQTL file. name_new- cohort name in the fixed eQTL file. add_to_file- whether or not add this cohort to the output (yes/no). NB! Pay attention that OverallZScore is still based on all the cohorts in original file! Order of the cohorts in this file is used to sort cohort ordering in the fixed file.')
                    
parser.add_argument('--orig_file', metavar = 'original eQTL mapping output file',
                    help = 'Original file form eQTLMapping output. Can be gzipped or unzipped.')
                    
parser.add_argument('--fixed_file', metavar = 'name of the fixed eQTL mapping output file',
                    help = 'Fixed file where cohort order and naming is updated. No file extensions needed and the file will be gzipped.')

parser.add_argument('--remove_cohort_specific_information',
                    help = 'If added, DatasetsWhereSNPProbePairIsAvailableAndPassesQC column will be replaced with the number of datasets; DatasetsZScores will be removed;	DatasetsNrSamples will be replaced with the sum of all the samples tested for this QTL', 
                    action = 'store_true')

args = parser.parse_args()

# Messages:
print 'Mapping file is:', args.mapping_file
print 'Original eQTL mapping output file:', args.orig_file
print 'Fixed eQTL mapping output file will be named:', args.fixed_file + '.txt.gz'
print 'Remove all cohort-specific information:', args.remove_cohort_specific_information

# Read in mapping file:
mapping = {}
keep = {}
cohort_order = []

with open(args.mapping_file, mode = 'r') as input_map:
  for line in input_map:
    line = line.strip('\n')
    (key, val, add_to_list) = line.split('\t')
    mapping[key] = val
    keep[key] = add_to_list
    cohort_order.append(key)


del mapping['name_old']
del keep['name_old']
cohort_order.pop(0)

# find a line which includes all cohorts:
# add .gz handle
if ".gz" in args.orig_file:
  print 'Input file is gzipped.'
  with gzip.open(args.orig_file, mode = 'r') as trans:
    line = trans.readline()
    line = trans.readline()
    
    while line:
      field = line.split('\t')
      field = field[11]
      
      line = trans.readline()
      
      if field.count(';-;') == 0:
        print "Assoc. with all the cohorts in the meta-analysis found!"
        break
      
    if field.count(';-;') > 0:
      print "Error! There are no eQTLs being tested in all the cohorts in the meta-analysis! Please check your analysis output!"
      quit()
else:
  print 'Input file is not gzipped.'
  with open(args.orig_file, mode = 'r') as trans:
    line = trans.readline()
    line = trans.readline()
    
    while line:
      field = line.split('\t')
      field = field[11]
      
      line = trans.readline()
      
      if field.count(';-;') == 0:
        print "Assoc. with all the cohorts in the meta-analysis found!"
        break
      
    if field.count(';-;') > 0:
      print "Error! There are no eQTLs being tested in all the cohorts in the meta-analysis! Please check your analysis output!"
      quit()

field = field.split(';')

#check if all the datasets are in the mapping file:
if len(set(cohort_order).intersection(field)) != len(set(field)):
  print "Error! There are fewer cohorts in the name mapping file than in the eQTL file!"
  quit()


# Remove elements not in the name mapping file:

keep = {k: v for k, v in keep.iteritems() if v == 'yes'}

# Make the index for sorting:
index_for_reorder = []

for i in range(0, len(field)):
    help_vect = field.index(cohort_order[i])
    index_for_reorder.append(help_vect)
# Sorting the full list:    
sorted_cohorts = []
for i in range(0, len(index_for_reorder)):
  sorted_cohorts.append(field[index_for_reorder[i]])

# Make index vector for keeping cohorts:
index_for_keeping = []
for i in range(0, len(sorted_cohorts)):
  if sorted_cohorts[i] in keep:
    help_vect = i
    index_for_keeping.append(help_vect)

################
# Testing bits #
################

# # test replacing with new name:
# for i in range(0, len(sorted_output)):
#   sorted_output[i] = mapping.get(sorted_output[i])
# 
# # test removing unwanted cohorts:
# filtered_cohorts = []
# for i in range(0, len(index_for_remove)):
#   filtered_cohorts.append(sorted_cohorts[index_for_remove[i]])

##################  
# Do the fixing: #
##################

with gzip.open(args.fixed_file + '.txt.gz', mode = 'w') as output:
  
  if ".gz" in args.orig_file:
    with gzip.open(args.orig_file, mode = 'r') as trans:
      
      line = trans.readline()
      # Change the column names if minimal output is needed:
      if args.remove_cohort_specific_information:
        fields = line.split('\t')
        fields[11] = 'NrOfCohortsTested'
        fields[12] = 'CohortZScores'
        fields[13] = 'SumNumberOfSamples'
        del fields[12]
        line = '\t'.join(fields)
        output.write(line)
      else:
        output.write(line)

      line = trans.readline()
      while line:
        # Split to fields:
        fields = line.split('\t')
        
        # Take out relevant fields (12, 13, 14):
        
        #################
        # Cohort names  #
        #################
        
        DatasetsWhereSNPProbePairIsAvailableAndPassesQC = fields[11]
        DatasetsWhereSNPProbePairIsAvailableAndPassesQC = DatasetsWhereSNPProbePairIsAvailableAndPassesQC.split(';')
        
        # if cohort-specific info removed then sum the sample sizes
        if args.remove_cohort_specific_information:
          DatasetNr = DatasetsWhereSNPProbePairIsAvailableAndPassesQC
          # If "-" then replace with 0 else replace with 1
          DatasetNr = [1 if x != '-' else x for x in DatasetNr]
          DatasetNr = [0 if x == '-' else x for x in DatasetNr]
          DatasetNr = [int(i) for i in DatasetNr]
          DatasetNr = sum(DatasetNr)
          fields[11] = str(DatasetNr)
          
        else:
          sorted_output = []
          
          for i in range(0, len(index_for_reorder)):
            sorted_output.append(DatasetsWhereSNPProbePairIsAvailableAndPassesQC[index_for_reorder[i]])
            
          # Update the naming of the cohorts
          for i in range(0, len(sorted_output)):
            if sorted_output[i] == '-':
              sorted_output[i] = '-'
            else:
              sorted_output[i] = mapping.get(sorted_output[i])
              
          
          # Removing unwanted cohorts:
          filtered_cohorts = []
          for i in range(0, len(index_for_keeping)):
            filtered_cohorts.append(sorted_output[index_for_keeping[i]])
            
          fields[11] = ';'.join(filtered_cohorts)
        
        #####################
        # Dataset Z-scores  #
        #####################
        
        DatasetsZScores = fields[12]
        DatasetsZScores = DatasetsZScores.split(';')
        
        sorted_output = []
        
        for i in range(0, len(index_for_reorder)):
          sorted_output.append(DatasetsZScores[index_for_reorder[i]])
          
        # Removing unwanted cohorts:
        filtered_cohorts = []
        for i in range(0, len(index_for_keeping)):
          filtered_cohorts.append(sorted_output[index_for_keeping[i]])
          
        fields[12] = ';'.join(filtered_cohorts)
        
        ###############
        # Dataset N's #
        ###############
        
        DatasetsNrSamples = fields[13]
        DatasetsNrSamples = DatasetsNrSamples.split(';')
        
        sorted_output = []
        for i in range(0, len(index_for_reorder)):
          sorted_output.append(DatasetsNrSamples[index_for_reorder[i]])
        
        # if cohort-specific info removed then sum the sample sizes
        if args.remove_cohort_specific_information:
          sorted_output = [0 if x == '-' else x for x in sorted_output]
          sorted_output = [int(i) for i in sorted_output]
          sorted_output = sum(sorted_output)
          fields[13] = str(sorted_output)
        else:
          # Removing unwanted cohorts:
          filtered_cohorts = []
          for i in range(0, len(index_for_keeping)):
            filtered_cohorts.append(sorted_output[index_for_keeping[i]])
          
          fields[13] = ';'.join(filtered_cohorts)

        if args.remove_cohort_specific_information:
          del fields[12]
          
        line_out = '\t'.join(fields)
        
        output.write(line_out)
        line = trans.readline()
        
  else:
    with open(args.orig_file, mode = 'r') as trans:
      
      line = trans.readline()
      # Change the column names if minimal output is needed:
      if args.remove_cohort_specific_information:
        fields = line.split('\t')
        fields[11] = 'NrOfCohortsTested'
        fields[12] = 'CohortZScores'
        fields[13] = 'SumNumberOfSamples'
        del fields[12]
        line = '\t'.join(fields)
        output.write(line)
      else:
        output.write(line)
        
      line = trans.readline()
      while line:
        # Split to fields:
        fields = line.split('\t')
        
        # Take out relevant fields (12, 13, 14):
        
        #################
        # Cohort names  #
        #################
        
        DatasetsWhereSNPProbePairIsAvailableAndPassesQC = fields[11]
        DatasetsWhereSNPProbePairIsAvailableAndPassesQC = DatasetsWhereSNPProbePairIsAvailableAndPassesQC.split(';')
        
        # if cohort-specific info removed then sum the sample sizes
        if args.remove_cohort_specific_information:
          DatasetNr = DatasetsWhereSNPProbePairIsAvailableAndPassesQC
          # If "-" then replace with 0 else replace with 1
          DatasetNr = [1 if x != '-' else x for x in DatasetNr]
          DatasetNr = [0 if x == '-' else x for x in DatasetNr]
          DatasetNr = [int(i) for i in DatasetNr]
          DatasetNr = sum(DatasetNr)
          fields[11] = str(DatasetNr)
          
        else:
          sorted_output = []
          
          for i in range(0, len(index_for_reorder)):
            sorted_output.append(DatasetsWhereSNPProbePairIsAvailableAndPassesQC[index_for_reorder[i]])
            
          # Update the naming of the cohorts
          for i in range(0, len(sorted_output)):
            if sorted_output[i] == '-':
              sorted_output[i] = '-'
            else:
              sorted_output[i] = mapping.get(sorted_output[i])
              
          
          # Removing unwanted cohorts:
          filtered_cohorts = []
          for i in range(0, len(index_for_keeping)):
            filtered_cohorts.append(sorted_output[index_for_keeping[i]])
            
          fields[11] = ';'.join(filtered_cohorts)
        
        #####################
        # Dataset Z-scores  #
        #####################
        
        DatasetsZScores = fields[12]
        DatasetsZScores = DatasetsZScores.split(';')
        
        sorted_output = []
        
        for i in range(0, len(index_for_reorder)):
          sorted_output.append(DatasetsZScores[index_for_reorder[i]])
          
        # Removing unwanted cohorts:
        filtered_cohorts = []
        for i in range(0, len(index_for_keeping)):
          filtered_cohorts.append(sorted_output[index_for_keeping[i]])
          
        fields[12] = ';'.join(filtered_cohorts)
        
      ###############
        # Dataset N's #
        ###############
        
        DatasetsNrSamples = fields[13]
        DatasetsNrSamples = DatasetsNrSamples.split(';')
        
        sorted_output = []
        for i in range(0, len(index_for_reorder)):
          sorted_output.append(DatasetsNrSamples[index_for_reorder[i]])
        
        # if cohort-specific info removed then sum the sample sizes
        if args.remove_cohort_specific_information:
          sorted_output = [0 if x == '-' else x for x in sorted_output]
          sorted_output = [int(i) for i in sorted_output]
          sorted_output = sum(sorted_output)
          fields[13] = str(sorted_output)
        else:
          # Removing unwanted cohorts:
          filtered_cohorts = []
          for i in range(0, len(index_for_keeping)):
            filtered_cohorts.append(sorted_output[index_for_keeping[i]])
          
          fields[13] = ';'.join(filtered_cohorts)

        if args.remove_cohort_specific_information:
          del fields[12]
          
        line_out = '\t'.join(fields)
        
        output.write(line_out)
        line = trans.readline()

print "Fixing has finished!"
