from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
import numpy as np
# create a new Firefox session
### WEBSCRAPE
caps = DesiredCapabilities().FIREFOX
caps["marionette"] = True
caps["pageLoadStrategy"] = "normal"  #  complete
driver = webdriver.Firefox(capabilities=caps)
driver.maximize_window()
page = 'http://crdd.osdd.net/raghava/cppsite/browse_sub1.php?token=Linear&col=5'
driver.get(page)
f1= open('0.txt','w')
f1.write(driver.page_source.encode('utf-8'))
f1.close()
for i in range(2,37):
    page = 'http://crdd.osdd.net/raghava/cppsite/browse_sub1.php?token=Linear&col=5&page='
    page = page + str(i)
    driver.get(page)
    f1= open('%s.txt'%i,'w')
    f1.write(driver.page_source.encode('utf-8'))
    f1.close()

### NOW READ TXT FILE ###
import pandas as pd
data = pd.DataFrame(columns=['text'])
import time
for i in range(2,37):
    f1=open(str(i)+'.txt','r')
    for line in f1:
        if 'a href="display_seq.php?details=' in line:
            line1 =  BeautifulSoup(line, "html5lib").get_text().encode('utf-8')
            line2 = line1.split('     ID')
            line2[0] = line2[0].split('IDID')[1]
            data = data.append(line2)
data = data.reset_index(drop=1)
def get_import_stuff(str):
    result = ['',]*13
    str1 = str.split('PEPTIDE SEQUENCE')
    result[0] = str1[0]
    change = ['PEPTIDE NAME','LENGTH','CHIRALITY','LINEAR/CYCLIC','SOURCE','CATEGORY',
              'N TERMINAL MODIFICATION','C TERMINAL MODIFICATION','CHEMICAL MODIFICATION','CARGO','PUBMED ID','\n']
    for i in range(0,len(change)):
        str1 = str1[1].split(change[i])
        result[i+1] = str1[0]
    return result

impt = ['ID','PEPTIDE SEQUENCE','PEPTIDE NAME','LENGTH','CHIRALITY','LINEAR/CYCLIC','SOURCE','CATEGORY',
        'N TERMINAL MODIFICATION','C TERMINAL MODIFICATION','CHEMICAL MODIFICATION','CARGO','PUBMED ID']
a=data[0].map(get_import_stuff)
for i in impt:
    data[i] = 0
for i in range(len(data)):
    data= data.set_value(i,impt,a[i])
    
data['seq'] = data['PEPTIDE SEQUENCE'].apply(lambda x : x.upper())
data.to_csv('peptide_CPP2.csv',index=0)
