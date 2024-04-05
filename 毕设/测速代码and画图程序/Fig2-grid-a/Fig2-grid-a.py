import matplotlib.pyplot as plt
##import numpy as np
#import matplotlib as mpl
import re
import matplotlib.gridspec as gridspec
#import matplotlib.ticker as mtick

#plt.rcParams['font.family'] = ['sans-serif']#
plt.rcParams['font.sans-serif'] = ['SimHei']

def getY(fp,yl,y_pat,sc):
    with open(fp,'r',encoding='utf-8')as f:
        for i in range(sc+1):
            yl.append(float(y_pat.search(f.readline()).group('s')))
            f.readline()

def getX(y,x):
    j=-1
    for item in y:
        j=j+1
        x.append(j*2)


scale=50
file_name1='./file1-a.txt'
file_name2='./file1-b.txt'

#fig=plt.figure()

######################### figure layout set #########
fig=plt.figure(figsize=[10,5])
#fig.subplots_adjust(hspace=0.3, wspace=0.3)

top1=0.96
height1=0.425

top2=top1-height1
height2=0.425

gs11 = gridspec.GridSpec(1, 1) #(col,row)
gs11.update(bottom=top1-height1,top=top1,left=0.08, right=0.96) #gs1.update(bottom=0.48,top=0.8)
ax11 = plt.subplot(gs11[0, 0])
#ax2 = plt.subplot(gs1[1, 0])

gs12 = gridspec.GridSpec(1, 1)
gs12.update(bottom=top2-height2,top=top2,left=0.08, right=0.96)
ax12 = plt.subplot(gs12[0, 0])
#ax5 = plt.subplot(gs2[1, 0])
##################################
#ax=fig.add_subplot(221)

#para_pat=re.compile('PMFS=(.{12})')  # . 表示任意字符，12表示12个
#para_pat=re.compile('PMFS=(.+?),.+') # + 表示一个或多个重复，？表示‘不贪婪’，即向前匹配
#para_pat=re.compile('nr=(?P<s>.+?),.+')
# 这三种都可以，前两种引用用group(1);最后一种，用？P来表示模式，可以用<>中的字符串引用此group,见上

para_pat=re.compile('cgr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
x=[]
getX(y,x)
ax11.plot(x,y,marker='o',markersize=4,ls='',c='b',label='Cgr')

para_pat=re.compile('cnr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax11.plot(x,y,marker='o',markersize=4,ls='',c='r',label='Cnr')

para_pat=re.compile('cnrgr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax11.plot(x,y,marker='o',markersize=4,ls='',c='purple',label='Cgr+nr')

ax11.text(0,190,'a', fontsize=18)

ax11.legend(bbox_to_anchor=(1,0.6), numpoints=1, fontsize=10)

ax11.set_ylim(-20, 220)

##################################
#ax=fig.add_subplot(223)

para_pat=re.compile('gr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax12.plot(x,y,c='b',linewidth=0.9,label='gr')

para_pat=re.compile('nr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax12.plot(x,y,c='r',linewidth=0.9,label='nr')

para_pat=re.compile('contr1=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax12.plot(x,y,c='g',linewidth=0.9,label='ctl')

para_pat=re.compile('nr=(?P<s>.+?),.+')
y=[]
getY(file_name1,y,para_pat,scale)
ax12.plot(x,y,c='r',linewidth=0.9)

ax12.legend(bbox_to_anchor=(1,0.55), numpoints=1, fontsize=10)

ax12.set_ylim(-100, 1600)

ax12.set_xlabel('Time ($10^4$ Monte Carlo Steps)')
'''
##################################
ax=fig.add_subplot(222)

#para_pat=re.compile('PMFS=(.{12})')  # . 表示任意字符，12表示12个
#para_pat=re.compile('PMFS=(.+?),.+') # + 表示一个或多个重复，？表示‘不贪婪’，即向前匹配
#para_pat=re.compile('nr=(?P<s>.+?),.+')
# 这三种都可以，前两种引用用group(1);最后一种，用？P来表示模式，可以用<>中的字符串引用此group,见上

para_pat=re.compile('cnr=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
x=[]
getX(y,x)
ax.plot(x,y,marker='o',markersize=4,ls='',c='r')

para_pat=re.compile('cgr=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
ax.plot(x,y,marker='o',markersize=4,ls='',c='b')

para_pat=re.compile('cnrgr=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
ax.plot(x,y,marker='o',markersize=4,ls='',c='purple')

ax.text(0,190,'b', fontsize=18)

ax.set_ylim(-20, 220)

##################################
ax=fig.add_subplot(224)

para_pat=re.compile('nr=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
ax.plot(x,y,c='r',linewidth=0.9)

para_pat=re.compile('contr1=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
ax.plot(x,y,c='g',linewidth=0.9)

para_pat=re.compile('gr=(?P<s>.+?),.+')
y=[]
getY(file_name2,y,para_pat,scale)
ax.plot(x,y,c='b',linewidth=0.9)

ax.set_ylim(-100, 1600)

ax.set_xlabel('Time ($10^4$ Monte Carlo Steps)')
'''
##################################
plt.savefig('Fig2-a.png',dpi=600)
plt.show()


