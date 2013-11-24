# -*- coding: utf-8 -*-

'''
Created on 27.09.2012

@author: Esir Pavel
'''

import numpy as np
import matplotlib.pyplot as pl
import csv

f = open('oscill.csv', "r")
rdr = csv.reader(f,delimiter=";")
t = []
V_mean = []
for l in rdr:
    t.append(l[0])
    V_mean.append(l[1])

pl.figure()
pl.plot(t, V_mean, 'r')
pl.title(u"Электроэнцефалограмма")
pl.ylabel(u'Средний мембранный потенциал, мВ')
pl.grid()
pl.xlabel(u"Время, мс")
pl.show()
