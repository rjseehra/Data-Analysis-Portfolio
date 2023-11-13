import time
import numpy as np
import webdrivermanager

from selenium import webdriver
from webdrivermanager import ChromeDriverManager
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options

#LOCATE CHROMEDRIVER RELATIVE TO SCRIPT FILE LOCATION
option = webdriver.ChromeOptions()
option.add_argument("--user-data-dir=C:/Path/To/User/Google/Chrome/User Data")

#SELECT SELENIUM OPTIONS (USER PROFILE, SANDBOX MODE, ETC.)
option.add_argument('--profile-directory=Profile 3') 
option.add_argument('--no-sandbox')
option.add_argument('-disable-dev-shm-usage')
option.binary_location=""

driver = webdriver.Chrome('chromedriver.exe', options=option)

#CALL URL
driver.get("https://login.myameego.com/")

#OBTAIN LOGIN ELEMENTS
driver.find_element(By.ID, "username")
driver.find_element(By.ID, "password")
driver.find_element(By.ID, "client-id")

#CREATE STORAGE ELEMENTS FOR LOGIN CREDENTIALS
letusername = driver.find_element(By.ID, "username")
letpassword = driver.find_element(By.ID, "password")
letclientid = driver.find_element(By.ID, "client-id")

#OBTAIN XPATH FOR LOGIN BUTTON
login = driver.find_element(By.XPATH, "//button[@type='submit']")

#STORE AND ENTER LOGIN CREDENTIALS
letusername.send_keys("username")
letpassword.send_keys("password")
letclientid.send_keys("client-id")

#CLICK SUBMIT / LOGIN AND LET LOAD 
login.click()
time.sleep(2)

shifts = []

#FIND SHIFT LIST BY COMMON CLASS VIA XPATH
element = driver.find_elements(By.XPATH, "//div[@class='col-sm-12']")

#TURN WEB ELEMENTS INTO READABLE TEXT
for i in element:
    shift = str(i.text).split('\n')
    shifts.append(shift)
del shifts[-1]

#CREATE NUMPY ARRAY TO STORE AND ORGANIZE SHIFT INFORMATION 
shiftarray = np.array(shifts)
for i in shiftarray:
    day = i[0]
    date = i[1]
    time = i[2]
  
#RETURN ORGANIZED SHIFT INFORMATION 
print(shiftarray)

#CLOSE CHROMEDRIVER
driver.close() 
