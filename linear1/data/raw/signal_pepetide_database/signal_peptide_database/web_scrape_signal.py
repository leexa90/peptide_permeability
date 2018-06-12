from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
import numpy as np
import time
import os
import pandas as pd


### WEBSCRAPE
##caps = DesiredCapabilities().FIREFOX
##caps["marionette"] = True
##caps["pageLoadStrategy"] = "normal"  #  complete
##driver = webdriver.Firefox()
##driver.maximize_window()
##page = 'http://www.signalpeptide.de/?m=searchspdb'
##driver.get(page)
##die

data_file = [x for x in  os.listdir('.') if x.endswith('.html')]
data = pd.DataFrame(columns = ['ID','entry name','protein name','organism','len','status','seq',])
for file in sorted(data_file):
    temp = ''
    f1 = open(file,'r')
    for line in f1:
        temp += line
    text = BeautifulSoup(temp, "html5lib").get_text()
    start = False
    result = []
    for line in text.split('\n'):
        if start:
            result += [line,]
        if 'Accession NumberEntry NameProtein NameOrganismLengthStatusSignal Sequence' in line:
            start = True
    for i in range(0,len(result)-5,1):
        if result[i] != '' and result[i+5] =='confirmed':
            newrow =pd.DataFrame([result[i:i+7]],columns = ['ID','entry name','protein name','organism','len','status','seq',])
            data = data.append( newrow)
    print file,len(data)

data.to_csv('Mammal_signal_peptide.csv',index=0)

