
# coding: utf-8

# In[1]:


# imports 
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


# Get xml file here: https://www.ncbi.nlm.nih.gov/biosample/?term=antibiogram%5bfilter%5d
xmlfile = '/home/ben/bhatt_local/meta-mustache/biosample_result.xml'
tree = ET.parse(xmlfile)
root = tree.getroot()


# In[3]:


# testing xml stuff... 
if False:
    for child in root[0:5]:
        print (child.tag, child.attrib)
    
    for i in root.iter('Table'):
        xmltable_to_df(i)
        for j in i.iter('Header'):
            for child in j:
                print(child.text)


# In[125]:


a=root[0][0]
for child in a:
    print(child.tag, list(child.attrib.values()))


# In[148]:


def parse_ncbi_xml(root):
    """return a tabel form parsing the xml file from NCBI
    Each row of the table is an enry about resistance to antibiotics
    for a particular organism
    """
    
    ii = 0
    bigdf = pd.DataFrame()
    for i in root:
        # basic sample information 
        sample_id = i.attrib['id']
        sample_title = i[1][0].text
        taxonomy_name = i[1][1].attrib['taxonomy_name']
        try:
            organism_name = i[1][1][0].text
        except IndexError:
            organism_name = 'NA'
        # check for sequencing data
        has_sra = False
        for child in i[0]:
            if 'SRA' in list(child.attrib.values()):
                has_sra =True
        try:
            comment_text = i[1][2].findall('Paragraph')[0].text
        except IndexError:
            comment_text = 'NA'
        if has_sra:
            # add this to each row
            add_extra = [sample_id, sample_title, taxonomy_name, organism_name, comment_text]
            add_extra_header = ['sample_id', 'sample_title', 'taxonomy_name', 'organism_name', 'comment_text']
            # get resistance table
            table = i[1][2].findall('Table')[0]
            table_title = table.attrib['class']
            table_header = [col.text for col in table.find('Header').getchildren()] + add_extra_header
            rows = [c for c in table.find('Body').getchildren()]
            table_text = []
            for r in rows:
                row_text = [cell.text for cell in r.getchildren()] + add_extra
                table_text.append(row_text)
            table_df = pd.DataFrame(table_text, columns=table_header)

            # standardize the columns in this 
            # we want the Antibiotic, Resistance phenotype, Measurement, Measurement sign, Measurement units
            # columns. If theres not that, put NA
            keep_columns = ['Antibiotic', 'Resistance phenotype', 'Measurement', 'Measurement sign', 'Measurement units'] + add_extra_header
            for k in keep_columns:
                if k not in table_df.columns:
                    table_df[k] = 'NA'
            table_df = table_df[keep_columns]
            # add other info to the table
    #         table_df['sample_id']= sample_id
    #         table_df['sample_title']= sample_title
    #         table_df['taxonomy_name']= taxonomy_name
    #         table_df['organism_name']= organism_name
            bigdf = bigdf.append(table_df)
        ii +=1
    return(bigdf)


# In[149]:


# get this big table
bigdf = parse_ncbi_xml(root)


# In[150]:


print(bigdf.shape)
# rename the columns
bigdf.columns = ['antibiotic', 'resistance_phenotype', 'measurement', 'measurement_sign',
                 'measurement_units', 'sample_id', 'sample_title', 'taxonomy_name', 'organism_name', 'comment_text']
        
bigdf.to_pickle('/home/ben/bhatt_local/meta-mustache/resistance_table.pickle')


# In[151]:


bigdf = pd.read_pickle('/home/ben/bhatt_local/meta-mustache/resistance_table.pickle')
print(bigdf.shape)


# In[152]:


antibiotic_freq = pd.crosstab(bigdf.antibiotic, 'antibiotic').sort_values('antibiotic', ascending=False)
antibiotic_freq[0:9]


# In[153]:


# is carbapanem in this list?
print('carbapenem' in antibiotic_freq)
# NOPE


# In[154]:


bigdf.groupby('taxonomy_name')['sample_id'].nunique().sort_values(ascending=False)[0:9]


# In[155]:


pd.crosstab([bigdf.antibiotic, bigdf.resistance_phenotype], bigdf.organism_name)


# In[156]:


ecolidf = bigdf.loc[bigdf.organism_name=='Escherichia coli']
pd.crosstab(ecolidf.resistance_phenotype, ecolidf.antibiotic, margins=True)


# In[157]:


ecoli_ciprodf = ecolidf.loc[ecolidf.antibiotic=='ciprofloxacin']


# In[158]:


pd.crosstab(ecoli_ciprodf.resistance_phenotype, ecoli_ciprodf.antibiotic)


# In[159]:


ecoli_ciprodf.measurement = pd.to_numeric(ecoli_ciprodf.measurement)


# In[182]:


ecoli_ciprodf = ecoli_ciprodf.sort_values('measurement')
susceptible_df = ecoli_ciprodf.loc[ecoli_ciprodf.resistance_phenotype=='susceptible']
susceptible_ids = susceptible_df['sample_id']
resistant_df = ecoli_ciprodf.loc[ecoli_ciprodf.resistance_phenotype=='resistant']
resistant_ids = resistant_df['sample_id']


# In[183]:


print(len(np.unique(susceptible_ids.values)))
print(len(np.unique(resistant_ids.values)))


# In[184]:


susceptible_ids.to_csv('/home/ben/bhatt_local/meta-mustache/ecoli_cipro_susceptible_ids.txt', index=False)
resistant_ids.to_csv('/home/ben/bhatt_local/meta-mustache/ecoli_cipro_resistant_ids.txt', index=False)
all_ids = pd.Series(np.append(susceptible_ids.values,resistant_ids.values))
all_ids.to_csv('/home/ben/bhatt_local/meta-mustache/ecoli_cipro_all_ids.txt', index=False)


# In[185]:


# save dataframe to go with this
# savedf = ecoli_ciprodf.iloc[0:50]
# savedf = savedf.append(ecoli_ciprodf.iloc[-50:])

savedf = susceptible_df.append(resistant_df)
savedf.to_csv('/home/ben/bhatt_local/meta-mustache/ecoli_cipro_sample_table.txt', index=False, sep='\t' )


# In[186]:


savedf.shape

