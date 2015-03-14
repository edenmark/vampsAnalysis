#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy
import matplotlib.pyplot as pyplot


""" Select the files in which you would like to extract data"""

filelist = ['noSeawater.txt']
metadata_list = []

""" Probability Distribution Functions and Classes"""

class DDist:
    def __init__(self, dictionary):
        if not (abs(sum(dictionary.values())-1) < 1e-6 and min(dictionary.values()) >= 0.0):
            raise Exception, "Probabilities must be nonnegative, and must sum to 1"
        self.d = dictionary

    def prob(self, elt):
        return self.d.get(elt,0)

    def support(self):
        #### POSSIBLE FUTURE BUG; SWITCHED TO >=         \/
        return [k for k in self.d.keys() if self.prob(k) >= 0]

    def draw(self):
        r = random.random()
        sum = 0.0
        for val in self.support():
            sum += self.prob(val)
            if r < sum:
                return val

    def __repr__(self):
        return "DDist(%s)" % repr(self.d)
    
    __str__ = __repr__

    def expectation(self,f):
        return sum(self.prob(i)*f(i) for i in self.support())
   
    def mean(self):
        return self.expectation(lambda x: x)

    def variance(self):
        m = self.mean()
        return self.expectation(lambda x: (x-m)**2)

    def stddev(self):
        return self.variance()**0.5

    def maxProbElt(self):
        return max(self.support(), key=self.prob)

    def project(self, mapFunc):
        result = {}
        for e in self.support():
            new = mapFunc(e)
            result[new] = result.get(new,0) + self.prob(e)
        return DDist(result)

    def condition(self, testFunc):
        newElements = [e for e in self.support() if testFunc(e)]
        z = sum(self.prob(e) for e in newElements)
        result = {}
        for e in newElements:
            result[e] = result.get(e,0) + self.prob(e)/z
        return DDist(result)

def bayesRule(PA, PBgA, b):
    return makeJointDistribution(PA, PBgA).condition(lambda x: x[1]==b).\
        project(lambda x: x[0])

def makeJointDistribution(PA, PBgA):
    d = {}
    for a in PA.support():
        for b in PBgA(a).support():
            d[(a, b)] = PA.prob(a) * PBgA(a).prob(b)
    return DDist(d)

def totalProbability(PA, PBgA):
    return makeJointDistribution(PA, PBgA).project(lambda x: x[1])



""""Creation of Dictionaries from extracted information"""

def make_dicts(filelist):
    """
    Creates and returns 5 dictionaries.
    Sample Dict: maps the number of each classification to sample
    Study Dict: maps total number of each classification to corresponding study
    Total Dict: maps total number of bacteria in each sample to each sample
    Total Study Dict: maps total number of bacteria to each study
    Total Bacteria Dict: how many of each bacteria there are in total
    """
    sample_dict = {}
    study_dict = {}
    total_dict = {}
    total_study_dict = {}
    total_bacteria_dict = {}

    for new_file in filelist:
        
        matrix = []
        new_file = open(new_file, 'r',0)
        for row in new_file:
            row = row.split('\t')
            row[len(row)-1]=row[(len(row)-1)].rstrip('\n')
            matrix.append(row)
            # we now have huge matrix representing the table
        first = True
        matrix = matrix[:(len(matrix)-6)]
        
        old_samples = []
        new_samples = []
        for sample in matrix[0]:
            '''
            Creates a dictionary to describe each sample and its bacteria distribution.
            In the form: {sample{bacteria1:number, bacteria2: number}}
            '''
            old_samples.append(sample.lstrip().rstrip('\n'))
            second = True
            if first == True:
                first = False
            else:
                old_sample = sample
                new_sample = sample.split('--')
                new_sample = new_sample[1]
                new_sample = new_sample.lstrip().rstrip('\n')
                new_samples.append(new_sample)
                sample_dict[new_sample] = {}
                col = matrix[0].index(old_sample)
                for each_row in matrix:
                    if second == True:
                        second = False
                    else:
                        row = matrix.index(each_row)
                        category = each_row[0]
                        sample_dict[new_sample][category] = int(matrix[row][col].rstrip('\n'))

        
        old_samples = old_samples[1:]
        n = 0
        for row in matrix:
            if matrix.index(row) == 0:
                pass
            else:
                name = row[0]
                total = 0
                for col in row[1:]:
                    col = col.rstrip('\n')
                    total += int(col)
                total_bacteria_dict[name] = total

        for sample in sample_dict:
            '''
            Creates dictionary of each sample and total amount of bacteria in each sample
            {sample:total_bacteria_num}
            '''
            total = 0
            for x in sample_dict[sample]:
                total+=sample_dict[sample][x]
            total_dict[sample] = total
        """
        for sample in old_samples:
            '''
            Creates a dictionary that concatenates the sample dictionary into just studies.
            Allows us to look at the study as a whole instead of each sample within the study.
            In the form:    {study:{bacteria1:num}}
            '''
            temp = sample.split('--')
            study = temp[0]
            sample = temp[1]
            # MUTABLE DICTIONARIES **********************************************************************
            if study not in study_dict:
                study_dict[study] = sample_dict[sample]
            else:
                for bacteria in sample_dict[sample]:
                    study_dict[study][bacteria] += sample_dict[sample][bacteria]
        
        for study in study_dict:
            '''
            Creates dictionary of each sample and total amount of bacteria in each study
            {study:total_bacteria_num}
            '''
            total = 0
            
            for x in study_dict[study]:
                total+=study_dict[study][x]
            total_study_dict[study] = total
        """
    
    
    return sample_dict, total_dict, total_bacteria_dict

dictionary = make_dicts(filelist)
sample_dict = dictionary[0]
#study_dict = dictionary[1]
total_dict = dictionary[1]
#total_study_dict = dictionary[3]
total_bacteria_dict = dictionary[2]
bact_list = []
for b in total_bacteria_dict:
    bact_list.append(b)
bact_list.sort()


"""
def create_large_file(s=sample_dict, t = total_bacteria_dict, b=bact_list):
     Reverses the dictionaries made into a text file
    s_list = []
    for q in sample_dict:
        s_list.append(q)
    s_list.sort()
    new_file = open('combined.txt', 'w', 0)
    firstline = 'taxa'
    n=0
    for sample in s_list:
        firstline += '\t' + str(sample)
    new_file.write(firstline + '\n')
    for bact in b:
        line = str(bact) + '\t'
        n +=1
        for samp in s_list:
            
            
            if bact not in s[samp]:
                line += str(0) + '\t'
            else:
                line += str(s[samp][bact]) + '\t'
        line += '\n'
        new_file.write(line)
    new_file.close()
"""
#create_large_file()
make_dicts(filelist)




def rip_metadata(metadata_list):
    """
    Creates and returns a dictionary of the metadata (species of coral) **SUBMITTED IN PROPER FORM** 
    Maps sample to coral species
    """
    metadata_dict = {}
    for filename in metadata_list:
        matrix = []
        new_file = open(filename, 'r', 0)
        for row in new_file:
        
            row = row.strip('\n')
            row = row.split('\t') 
            matrix.append(row)

        for row in matrix:
            if row[1] == 'sample_id':
                x= matrix.index(row)
        n = 4
        for sample in matrix[0][4:]:
            metadata_dict[sample] = matrix[x][n]
            n+=1
    return metadata_dict

metadata_dict = rip_metadata(['cmp1.txt','cmp2.txt','cmp3.txt','cmp4.txt','cmp5.txt','cmp6.txt','cmp7.txt'])


def removeNonAscii(s):
    return "".join(i for i in s if ((ord(i)<128) and ord(i)!=0))
                   
def getWeirdAIndex(s):
    for i in s:
        if ord(i) == 194:
            return s.index(i)
    return 0          

def create_species_dict(m = metadata_dict, s = sample_dict):
    """
    Of all the CORAL species in all of the studies, this will map each CORAL SPECIES to the number of each
    BACTERIAL SPECIES within each coral. 
    """
    species_dict = {}
    coraltosample = {}
    
    for each_sample in sample_dict:
        
        each_sample = each_sample.lstrip().rstrip('\n')
        
        if each_sample == 'WP':
            
            coral_species = m['MannWP']
        elif each_sample == 'ApDMu3':
            coral_species = m['ApMuD2']
        elif each_sample == 'ApDMu2':
            coral_species = m['ApMuD2']
        elif each_sample == 'ApDMu1':
            coral_species = m['ApMuD1']
        elif each_sample == 'ApNMu5':
            coral_species = m['ApMuN5']
        elif each_sample == 'ApNMu4':
            coral_species = m['ApMuN4']
        elif each_sample == 'ApNMu6':
            coral_species = m['ApMuN6']
        elif each_sample == 'SW17_CG2':
            coral_species = m['SW17CG2']
        elif each_sample == 'SW11_AQ3':
            coral_species = m['SW_11AQ3']
        else:
            coral_species = m[each_sample]
        coral_species = coral_species.lstrip().rstrip()

        index = getWeirdAIndex(coral_species)
        coral_species = removeNonAscii(coral_species)
        coral_species = (coral_species[:index] +" " + coral_species[index:]).lstrip()
        if coral_species not in coraltosample:
            coraltosample[coral_species] = [each_sample]
        else:
            coraltosample[coral_species].append(each_sample)

    # sample dict is messed up here in ApD1. dont know why
    k = 0
    
    for coral in coraltosample:
    
        species_dict[coral] = {}
        for sample in coraltosample[coral]:
            for bacteria in sample_dict[sample]:
                #print 'sample:', sample,'\t', ' bacteria:',  bacteria, '\t', ' num:', sample_dict[sample][bacteria]
                if bacteria not in species_dict[coral]:
                    species_dict[coral][bacteria] = sample_dict[sample][bacteria]
                else:
                    species_dict[coral][bacteria] += sample_dict[sample][bacteria]           
    return species_dict

species_dict = create_species_dict()



def get_percentage_sample(bacteria,sample,total=total_dict,d=sample_dict):
    if type(total) != int:
        total = total_dict[sample]
    num = d[sample][bacteria]
    return float(num)/total

def where_is_bacteria(bacteria, s = sample_dict, m = metadata_dict, p = species_dict, t = total_bacteria_dict):
    """
    Given a specified bacteria, returns a distribution of the bacteria across ALL CORAL SPECIES 
    """
    
    dist_dict = {}
    temp_d = {}
    for coral_species in p:
        
        
        for each_bacteria in p[coral_species]:
            if each_bacteria not in temp_d:
                temp_d[each_bacteria] = p[coral_species][each_bacteria]
            else:
                temp_d[each_bacteria] += p[coral_species][each_bacteria]
    
    total = total_bacteria_dict[bacteria]
    
    for coral in p:
        for bacteria_species in p[coral]:
            if bacteria_species == bacteria:
                num = p[coral][bacteria_species]
                #print 'coral: ', coral, ' bacteria: ', bacteria, ' num: ', num
                dist_dict[coral] = float(num)/total
    return DDist(dist_dict), bacteria

    
def whats_in_coral(coral, t = total_dict, m = metadata_dict, s = species_dict):
    """
    Given a specified coral species, returns the distribution of bacteria
    in that coral ACROSS ALL STUDIES
    """
    # find total num bact in coral
    # for each bact find num bact in coral
    # coral_num_dict will show how many total microbes in each coral species
    dist_dict = {}
    coral_num_dict = {}
    #print t
    theSum = 0
    for sample in t:
        
        if sample == 'WP':
            coral_species = m['MannWP']
        elif sample == 'ApDMu3':
            coral_species = m['ApMuD2']
        elif sample == 'ApDMu2':
            coral_species = m['ApMuD2']
        elif sample == 'ApDMu1':
            coral_species = m['ApMuD1']
        elif sample == 'ApNMu5':
            coral_species = m['ApMuN5']
        elif sample == 'ApNMu4':
            coral_species = m['ApMuN4']
        elif sample == 'ApNMu6':
            coral_species = m['ApMuN6']
        elif sample == 'SW17_CG2':
            coral_species = m['SW17CG2']
        elif sample == 'SW11_AQ3':
            coral_species = m['SW_11AQ3']
        else:
            coral_species = m[sample]

        # Getting NonAscii out of coral species names
        index = getWeirdAIndex(coral_species)
        coral_species = removeNonAscii(coral_species)
        coral_species = (coral_species[:index] +" " + coral_species[index:]).lstrip()

        if coral_species not in coral_num_dict:
            coral_num_dict[coral_species] = t[sample]
        else:
            coral_num_dict[coral_species] += t[sample]
        # coral num dict is made and describes total numbers of coral species
        
    total = coral_num_dict[coral]
    #print 'Total: ', total

    #print coral
    theSum = 0
    for bacteria in s[coral]:
        num = s[coral][bacteria]
        #print num, " " ,bacteria, " bacteria  in ", coral
        theSum += num
        #print 'theSum: ', theSum, num
        dist_dict[bacteria] = float(num)/total

    return DDist(dist_dict)


#bact_dist = where_is_bacteria('Bacteria;Acidobacteria')


def all_distributions(t=total_bacteria_dict, s=species_dict):
    """
    Returns a list of distributions of every bacteria
    """
    bact_dist_list = []
    coral_dist_list = []
    for bacteria in t:
        # create a distribution for each
        distribution = where_is_bacteria(bacteria)
        bact_dist_list.append(distribution)
    
    for coral_species in s:
        new_dist = whats_in_coral(coral_species)
        coral_dist_list.append(new_dist)
    
    return bact_dist_list , coral_dist_list



def compare_distributions(dist_list, t=total_bacteria_dict):
    """
    We now want to compare each one of our distributions to see if there is any
    similarities or differences in each bacterias distribution over coral species
    """

    # FIRST WE WILL TRY TO COMPARE EVERY SET OF BACTERIA TO SEARCH FOR SIMILARITIES
    # N! comparisions

    def similar(dist1, dist2, threshold):
        # returns True if dist1 and dist2 are 'similar' ; else False
        # makes the assumption that dist1 and dist 2 are distinctly different
        # see KL divergence
        total_dif = 0
        for coral in dist1.support():
            single_dif = abs(dist1.prob(coral)-dist2.prob(coral))
            total_dif += single_dif
        if total_dif < threshold:
            return True
        return False

    def inverse_correlation(dist1, dist2,threshold):
        # returns True if dist1 and dist2 are 'different' in that their 
        # presence inversely correlates with the presence of the other
        # makes the assumption that dist1 and dist 2 are distinctly different
        total_dif = 0
        for coral in dist1.support():
            single_dif = abs(dist1.prob(coral)-(1 - dist2.prob(coral)))
            total_dif += single_dif
        if total_dif < threshold:
            return True
        return False

    def uniform(dist1):
        # compares our distribution in question to a uniform distribution 
        num_coral = len(dist1.support())
        prob = 1.0/num_coral
        uni_dist = {}
        for coral in dist1.support():
            uni_dist[coral] = prob
        uniform_distribution = DDist(uni_dist)
        if similar(dist1,uniform_distribution,threshold):
            return True
        return False

    def peak(dist1, peakThreshold):
        # searches through our distribution to see if any single coral dominates bacterias distribution
        # returns a list of tuples specifying (coral, percent of bacteria distribution in specific coral)
        return_list = []
        for coral in dist1.support():
            if dist1.prob(coral) > peakThreshold:
                return_list.append((coral, dist1.prob(coral)))
        return return_list

    def visual(dist):
        
        pass

       
    threshold = len(t) * .025   ######### ARBITRARY THRESHOLD ##############
    already_seen = []           # this will be a list of tuples of previously seen dist pairs
    similar_list = []
    inverse_list = []
    uniform_list = []
    peak_list = []
    
    for tup1 in dist_list:
        dist = tup1[0]
        bact = tup1[1]
        for tup2 in dist_list:
            other_dist = tup2[0]
            other_bact = tup2[1]
            if dist == other_dist:
                pass
            else:
                if (bact, other_bact) not in already_seen:
                    if similar(dist,other_dist, threshold):                                       
                        similar_list.append((bact,other_bact))
                    if inverse_correlation(dist,other_dist,threshold):
                        inverse_list.append((bact,other_bact))
                    if uniform(dist):
                        uniform_list.append(bact)
                    if peak(dist,.4) != []:
                        peak_list.append((bact,peak(dist,.4)))
                    already_seen.append((bact,other_bact))
                    already_seen.append((other_bact,bact))

    ### IMPLEMENT CHECKS TO MAKE SURE NO 2 PAIRS ARE IN SIMILAR AND INVERSE

    for tup in similar_list:
        if tup in inverse_list:
            print 'WOOOOOAAHHHHH MAJOR ERROORRR'
    print 'phew made it'
    
     
    
    return similar_list, inverse_list, uniform_list, peak_list

dist_list = all_distributions()
compare_distributions(dist_list[0])


