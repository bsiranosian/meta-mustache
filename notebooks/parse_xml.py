
# coding: utf-8

# In[341]:


import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


tree = ET.parse('/home/ben/Downloads/biosample_result.xml')


# In[45]:


root = tree.getroot()


# In[9]:


for child in root[0:10]:
    print (child.tag, child.attrib)


# In[13]:


for child in root[0]:
    print (child.tag, child.attrib)


# In[25]:


root[0][1][1].attrib


# In[19]:


a.text


# In[26]:


get_ipython().run_line_magic('pinfo', 'root.get')


# In[37]:


for i in root.iter('Organism'):
    print(i.get('taxonomy_name'))


# In[54]:


for i in root.iter('Table'):
    xmltable_to_df(i)
    for j in i.iter('Header'):
        for child in j:
            print(child.text)


# In[153]:


a = root[42]
for child in a[1][2]:
    print (child.tag, child.attrib)


# In[154]:


atable = a[1][2][0]


# In[155]:


table_header = [c.text for c in atable.find('Header').getchildren()]
print(table_header)


# In[156]:


rows = [c for c in atable.find('Body').getchildren()]


# In[157]:


table_text = []
for r in rows:
    row_text = [cell.text for cell in r.getchildren()]
    table_text.append(row_text)


# In[158]:


table_text[0:3]


# In[159]:


df = pd.DataFrame(table_text, columns=table_header)


# In[160]:


df


# In[253]:


ii = 0
sample_dict = {}
for i in root:
    sample_id = i.attrib['id']
    sample_title = i[1][0].text
    taxonomy_name = i[1][1].attrib['taxonomy_name']
    try:
        organism_name = i[1][1][0].text
    except IndexError:
        organism_name = 'NA'
    table = i[1][2].findall('Table')[0]
#     if len(table) != 1:
#         print(ii)
    table_title = table.attrib['class']
#     print(str(ii))
#     print(str(ii) + ' ' + sample_id + ' ' + sample_title + ' ' + table_title)
    
    table_header = [col.text for col in table.find('Header').getchildren()]
    rows = [c for c in table.find('Body').getchildren()]
    table_text = []
    for r in rows:
        row_text = [cell.text for cell in r.getchildren()]
        table_text.append(row_text)
#     print(table_text)
#     print(table_header)
    
    table_df = pd.DataFrame(table_text, columns=table_header)
    
    sample_dict[sample_id] = [sample_id, sample_title, table_df, taxonomy_name, organism_name]
#     print (str(ii) + ' ' + table.tag)
    
#     if (i[1][2][0].tag) != 'Table':
#         print(str(ii) + '  ' + i[1][0].text)
    ii +=1


# In[190]:


d  =list(sample_dict.values())[0][2]
# any( =='amikacin')


# In[192]:


any(d.any =='amikacin')


# In[196]:


any(d=='amikacin')


# In[199]:


list(d.columns)


# In[214]:


# column headers of data table
col_headers = [list(i[2].columns) for i in list(sample_dict.values())]
col_headers_unlist = pd.Series([item for sublist in col_headers for item in sublist])
pd.crosstab(col_headers_unlist, columns='freq')
# all have the Antibiotic column


# In[235]:


antibiotics = [list(i[2]['Antibiotic']) for i in list(sample_dict.values())]
antibiotics_unlist = pd.Series([item for sublist in antibiotics for item in sublist])
antibiotics_table = pd.crosstab(antibiotics_unlist, columns='freq').sort_values('freq', ascending=False)
# top 5 occurances
antibiotics_table[0:5]


# In[236]:


# is carbapanem in this list?
print('carbapenem' in antibiotics_table)
# NOPE


# In[257]:


# Which taxonomies are in this list
tax_names = pd.Series([i[3] for i in list(sample_dict.values())])
tax_table = pd.crosstab(tax_names, columns='freq').sort_values('freq', ascending=False)
# top occurances
tax_table


# In[256]:


# organism name might be less restrictive
org_names = pd.Series([i[4] for i in list(sample_dict.values())])
org_table = pd.crosstab(org_names, columns='freq').sort_values('freq', ascending=False)
# top occurances
org_table
# not really. and some organisms are missing
# doesn't capture all the e coli strains


# In[309]:


# for e coli, get all antibiotic tables
my_organism = 'Escherichia coli'
ecoli_antibiotics = [list(i[2]['Antibiotic']) for i in list(sample_dict.values()) if i[3]==my_organism]
ecoli_antibiotics_unlist = pd.Series([item for sublist in ecoli_antibiotics for item in sublist])


# In[271]:


ecoli_antibiotics_table = pd.crosstab(ecoli_antibiotics_unlist, columns='freq').sort_values('freq', ascending=False)
ecoli_antibiotics_table[0:9]


# In[343]:


# for ecoli that have been tested with some antibiotic, what is the distribution of values?
ab = 'ciprofloxacin'
r= test_df.loc[test_df.Antibiotic =='ciprofloxacin', 'Resistance phenotype'].values[0]
resistance_bool_temp = [i[2].loc[i[2].Antibiotic == ab, 'Resistance phenotype'].values for i in 
                     list(sample_dict.values()) if i[3]=='Escherichia coli']
resistance_bool = [r[0] if len(r)>0 else 'NA' for r in resistance_bool_temp]
resistance_values_temp = [i[2].loc[i[2].Antibiotic == ab, 'Measurement'].values for i in 
                     list(sample_dict.values()) if i[3]=='Escherichia coli']
resistance_values = [float(r[0]) if len(r)>0 else np.nan for r in resistance_values_temp]
resistance_sign_temp = [i[2].loc[i[2].Antibiotic == ab, 'Measurement sign'].values for i in 
                     list(sample_dict.values()) if i[3]=='Escherichia coli']
resistance_sign = [r[0] if len(r)>0 else 'NA' for r in resistance_sign_temp]
resistance_unit_temp = [i[2].loc[i[2].Antibiotic == ab, 'Measurement units'].values for i in 
                     list(sample_dict.values()) if i[3]=='Escherichia coli']
resistance_unit = [r[0] if len(r)>0 else 'NA' for r in resistance_unit_temp]
resistance_text = [i1 +' '+str(i2)+' '+i3 for i1,i2,i3 in zip(resistance_sign, resistance_values, resistance_unit)]
ecoli_ids = [i[0] for i in list(sample_dict.values()) if i[3]=='Escherichia coli']
resistance_table = pd.DataFrame(data={'organism': my_organism, 'id': ecoli_ids, ab: resistance_bool,
                                      'measurement': resistance_text, 'measurement_value': resistance_values})


# In[371]:


resistance_table = resistance_table.sort_values(by='measurement_value', ascending=True, na_position='last')
resistance_table = resistance_table.loc[~np.isnan(resistance_table.measurement_value)]
resistance_table
# take first and last 5 ids
nrow = resistance_table.shape[0]
susceptible_ids = resistance_table.iloc[0:5]['id'].values
resistant_ids = resistance_table.iloc[nrow-5:]['id'].values


# In[353]:


rv = np.array(resistance_values)
plt.hist(rv[~np.isnan(rv)]); plt.show()


# In[373]:


print(resistant_ids)
print(susceptible_ids)


# In[1]:


1

