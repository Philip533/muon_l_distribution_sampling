import requests
from bs4 import BeautifulSoup

res1 = requests.get('https://xdb.lbl.gov/Section1/Table_1-1a.htm')
res2 = requests.get('https://xdb.lbl.gov/Section1/Table_1-1b.htm')
res3 = requests.get('https://xdb.lbl.gov/Section1/Table_1-1c.htm')
soup1 = BeautifulSoup(res1.content, 'html.parser')
soup2 = BeautifulSoup(res2.content, 'html.parser')
soup3 = BeautifulSoup(res3.content, 'html.parser')

tables1 = [
    [
        [td.get_text(strip=True) for td in tr.find_all('td')] 
        for tr in table.find_all('tr')
    ] 
    for table in soup1.find_all('table')
]
tables2 = [
    [
        [td.get_text(strip=True) for td in tr.find_all('td')] 
        for tr in table.find_all('tr')
    ] 
    for table in soup2.find_all('table')
]

tables3 = [
    [
        [td.get_text(strip=True) for td in tr.find_all('td')] 
        for tr in table.find_all('tr')
    ] 
    for table in soup3.find_all('table')
]

new_table1a = []
new_table1b = []
new_table2 = []
new_table3 = []
for i in range(1,len(tables1[0])):
    new_table1row = []
    for v in tables1[0][i]:
        if (v != ""):
            new_table1row.append(v)
    new_table1a.append(new_table1row)

for i in range(1,len(tables1[2])):
    new_table1row = []
    for v in tables1[2][i]:
        if (v != ""):
            new_table1row.append(v)
    new_table1b.append(new_table1row)

for i in range(1,len(tables2[0])):
    new_table2row = []
    for v in tables2[0][i]:
        if (v != ""):
            new_table2row.append(v)
    new_table2.append(new_table2row)
for i in range(1,len(tables3[0])):
    new_table3row = []
    for v in tables3[0][i]:
        if (v != ""):
            new_table3row.append(v)
    new_table3.append(new_table3row)

f = open("binding.dat", "a")
for i in range(len(new_table1a)):
    s = " ".join(map(str,new_table1a[i]))
    f.write(str(len(new_table1a[i][1:]))+" "+s+"\n")
for i in range(len(new_table1b)):
    s = " ".join(map(str,new_table1b[i]))
    f.write(str(len(new_table1b[i][1:]))+" "+s+"\n")
print("")
for i in range(len(new_table2)):
    s = " ".join(map(str,new_table2[i]))
    f.write(str(len(new_table2[i][1:]))+" "+s+"\n")
# print("")
for i in range(len(new_table3)):
    s = " ".join(map(str,new_table3[i]))
    f.write(str(len(new_table3[i][1:]))+" "+s+"\n")

