#!/usr/bin/python3

import os
import re
import tempfile
import lxml.etree
import pandas as pd
import numpy as np

# This script inspired by http://schoolofdata.org/2013/06/18/get-started-with-scraping-extracting-simple-tables-from-pdf-documents/

# source https://bitbucket.org/ScraperWiki/scraperwiki/raw/d7d079a44d8c1eccf07679d0ddb0de0b203178f1/scraperlibs/python/scraperwiki/utils.py
def pdftoetree(pdforig, options=''):
    """converts pdf file to xml file"""
    xmlin = tempfile.NamedTemporaryFile(mode='r', suffix='.xml')
    tmpxml = xmlin.name # "temph.xml"
    cmd = '/usr/bin/pdftohtml -xml -nodrm -zoom 1.5 -enc UTF-8 -noframes %s "%s" "%s"' % (options, pdforig, os.path.splitext(tmpxml)[0])
    cmd = cmd + " >/dev/null 2>&1" # can't turn off output, so throw away even stderr yeuch
    os.system(cmd)
    tree = lxml.etree.parse(xmlin)
    xmlin.close()
    return tree

def is_in_yrange(pos, ymin, ymax, trash):

    def lessthan(a, b):

        if a[0] < b[0]:
            return True
        elif a[0] == b[0] and a[1] < b[1]:
            return True
        else:
            return False

    def lessorequal(a, b):

        if a == b or lessthan(a, b):
            return True
        else:
            return False

    if lessorequal(ymin, pos) and lessthan(pos, ymax) and pos not in trash:
         return True
    else:
        return False

def extract_vals(pages, lims):
    """ Pulls out strings in certain region of page """

    data = list()
    for page in pages:
        pnum = int(page.attrib['number'])
        for el in page:
            if el.tag == 'text':
                xpos = int(el.attrib['left']) + int(el.attrib['width'])
                ypos = (pnum, int(el.attrib['top']))
                if is_in_yrange(ypos, **lims):
                    foo = lxml.etree.tostring(el,  method='text', encoding='unicode',
                                              with_tail=False)
                    foo = str.split(foo)
                    for i, val in enumerate(foo):
                        record = [ypos, xpos, i, val]
                        data.append(record)
    return data

def extract_limits(pages):
    """ Figures out region of pages containing table """

    lims = {}
    trash = []
    for page in pages:
        pnum = int(page.attrib['number'])
        for el in page:
            if el.tag == 'text':
                xpos = int(el.attrib['left']) + int(el.attrib['width'])
                ypos = int(el.attrib['top'])
                foo = lxml.etree.tostring(el, method='text', encoding='unicode', with_tail=False)
                pos = (pnum, ypos)
                if ypos > 870:
                    trash.append(pos)
                elif foo.find('Week') != -1 and 'ymin' not in lims:
                    lims['ymin'] = pos
                elif foo.find('Subtotal') != -1 and pos not in lims.values():
                    trash.append(pos)
                elif foo.find('TOTAL') != -1 and 'ymax' not in lims:
                    lims['ymax'] = pos
    lims['trash'] = set(trash)
    return lims

root_dir = '/home/docker/data/'
name = 'www.aasv.org_pedv_PEDV_weekly_report_140108.pdf'

pdfpath = root_dir + name

tree = pdftoetree(pdfpath)
pages = list(tree.getroot())

## Convert age-structured table

tabpages = [pages[1]]
lims = extract_limits(tabpages)
data = extract_vals(tabpages, lims)

df = pd.DataFrame(data, columns=['ypos', 'xpos', 'rank', 'val'])
dfs = df.sort(['ypos', 'xpos', 'rank'])

grouped = dfs.groupby('ypos', as_index=False)

dl = list()
for name,group in grouped:
    dl.append(' '.join(group.val).split())
header = dl[0]
assert header == ['Week', 'of', 'Submissions', 'Testing', 'Positive', 'for',
                  'PEDv', 'Suckling', 'Nursery', 'Grower', '/', 'Finisher',
                  'Sow', '/', 'Boar', 'Unk']
header = ['week','totalNumberSwineAccessions','Suckling','Nursery','Grower/Finisher','Sow/Boar','Unk']

dl = dl[1:]
df = pd.DataFrame(dl, columns=header)

df.to_csv('PEDvweeklyreport-age-ts-01-08-14.csv', index=False)

## Convert state-structured table

tabpages = pages[2:4]
lims = extract_limits(tabpages)
data = extract_vals(tabpages, lims)

df = pd.DataFrame(data, columns=['ypos', 'xpos', 'rank', 'val'])
dfs = df.sort(['ypos', 'xpos', 'rank'])

grouped = dfs.groupby('ypos', as_index=False)

dl = list()
for name,group in grouped:
    dl.append(' '.join(group.val).split())

header = dl[0]
assert header == ['Week', 'of', 'Positive', 'for', 'PEDv', 'CA', 'CO', 'IA',
                  'IL', 'IN', 'KS', 'KY', 'MD', 'MI', 'MN', 'MO', 'NC', 'NE',
                  'NY', 'OH', 'OK', 'PA', 'SD', 'TN', 'TX', 'WI', 'WY', 'Unk']
header = ' '.join(header)
header = re.sub('Week of', 'week', header)
header = re.sub('Positive for PEDv', 'totalNumberSwineAccessions', header)
header = header.split()

dl = dl[1:]

flag = False
dlTop = list()
dlBot = list()
for i in dl:
    if i[0] == '6/16/2013':
        flag = True
    if flag:
        dlBot.append(i)
    else:
        dlTop.append(i)

dfBot = pd.DataFrame(dlBot, columns=header)

botOnly = set(['MD', 'NE', 'CA', 'WY'])
headerTop = [lab for lab in header if lab not in botOnly]
dfTop = pd.DataFrame(dlTop, columns=headerTop)
df = pd.merge(dfTop, dfBot, how='outer')
df.to_csv('PEDvweeklyreport-state-ts-01-08-14.csv', index=False)

## Convert state-age-class cross-table

tabpages = [pages[4]]
lims = extract_limits(tabpages)
data = extract_vals(tabpages, lims)

df = pd.DataFrame(data, columns=['ypos', 'xpos', 'rank', 'val'])
dfs = df.sort(['ypos', 'xpos', 'rank'])

d1 = dfs[dfs.ypos >=(5, 150)]
state_col = d1[d1.xpos <= 180].val
state_col.name = 'State'

d2 = d1[d1.xpos >=600]
d2['col'] = np.nan

for i,xp in enumerate(d2.xpos):
    if xp <= 681:
        d2.iloc[i, 4] = 1
    elif xp > 681 and xp <= 747:
        d2.iloc[i, 4] = 2
    elif xp > 747 and xp <= 870:
        d2.iloc[i, 4] = 3
    elif xp > 870 and xp <= 959:
        d2.iloc[i, 4] = 4
    elif xp > 959 and xp <= 1059:
        d2.iloc[i, 4] = 5

d3 = d2[d2['rank'] == 0]

ind = list(d3.val == 'Grower/').index(True)
d3.iloc[ind, 3] = 'Grower/Finisher'
pv = d3.pivot(index='ypos', columns='col', values='val')

colnames = pv.iloc[0,:]
colnames = colnames.to_dict().values()

pv.columns = colnames
pv = pv[1:]
pv.index = state_col[1:]

pv[(np.isnan(pv.astype('float64')))] = 0
pv = pv.astype('int')

pv.to_csv('PEDvweeklyreport-state-age-cummulative-01-08-14.csv')


#print lxml.etree.tostring(pages[1], pretty_print=True)
