'''
author pgrbnv
http://mathprofi.ru/funkcia_raspredeleniya_dsv.html
http://mathprofi.ru/nepreryvnaya_sluchaynaya_velichina.html
'''
from math import *
import pylab 
from matplotlib import mlab
import sys
from multiprocessing import Process
import time
import warnings
from sympy import *
import sys
import string

warnings.filterwarnings('ignore')
__metaclass__ = type

class MethodError(Exception): pass

class plot_draw:
    def __init__(self,funct,xmin,xmax,dx,nx = False):
        self.xes = list(mlab.frange(xmin,xmax,dx))
        if nx != False:
            for m_y in nx:
                if m_y not in self.xes:
                    self.xes.append(m_y)
        self.yes = [funct(x) for x in self.xes]
    
    def show(self,style = ".g"):
        try :
            pylab.grid()
            pylab.plot(self.xes,self.yes,style)
            pylab.show()
        except AttributeError, e:
            print e
            print "fill the data"
        
class constructor(plot_draw) :
    def __init__(self):
        pass
        
    def __str__(self):
        return str(self.table)+"\nSee help()"
    
    @classmethod
    def decor(smth,a):
        def my_dec(func):
            def wrap(s,f=True):
                if f and type(f) == type(True):
                    print a[0],
                elif not f and type(f) == type(True):
                    print a[1],
                else:
                    print a[0],
                    
                print func(s,f)
            wrap._original = func
            return wrap
        return my_dec
    

    def graph(self):
        proc = []
        if self.flag :
            proc.append(Process(target=self.__dis_graph_proc,name='GraphDRV', args=()))
        else:
            proc.append(Process(target=self.__con_graph_proc,name='GraphCRV0', args=(True,)))
            proc.append(Process(target=self.__con_graph_proc,name='GraphCRV1', args=(False,)))
            
        for u in proc:
            u.start()
        time.sleep(0.1)
        
        
    def __dis_graph_proc(self):
        d = self.sortedvalues()
        ave = sum([d[f+1]-d[f] for f in range(len(d)-1)])/(len(d)-1)
        xm = (d[0]-ave,d[-1]+ave)
        f = self.function_top(d)
        super(constructor,self).__init__(f,xm[0],xm[1],ave/400.0,d)
        self.show()
    
    def __con_graph_proc(self,f):
        k = "-g"
        if not f:
            self.table = self.den._original(self)
            k = ".g"
        d = self.sortedvalues()
        ave = sum([(abs(float(u[0]))+abs(float(u[1])))/2 for u in self.table if u[0]!= -oo and u[1]!= oo])/(len(self.table.keys())-2)
        xm = (float(d[0][1]-ave),float(d[-1][0]+ave))
        f = self.ice_top(d)
        super(constructor,self).__init__(f,xm[0],xm[1],ave/400.0)
        self.show(k)
    
    def function_top(self,sx):
        tb = self.table
        def func(x):
            y = 0
            if x > sx[-1] :
                return 1
            if x < sx[0]:
                return 0 
            for u in range(len(sx)-1):
                if x > sx[u] and x <= sx[u+1]:
                    for q in sx:
                        if q == sx[u+1]:
                            break
                        y += tb[q]
            return y
        return func
    
    def ice_top(self,sx):
        inf = self.table
        def funct(x):
            for u in sx:
                if x > u[0] and x <= u[1]:
                    return inf[u][1].subs(symbols(inf[u][0]),x)
        return funct
        
    def sortedvalues(self):
        return sorted([s for s in self.table])
    
    def error_exeption(self,s,e):
        print "uncountable now"
        print "\n"
        raise MethodError(s+'('+str(e)+")")


class discrete(constructor):
    
    def __init__(self,dep):
        if type(dep) == type(dict()):
            self.table = dep
        else:
            print "error in initiation, dict expected"
    
    def help(self):
        print "Discrete Random Value"
        try :
            for i in self.table:
                print str(i)+" - "+str(self.table[i])
            print "dispersion -> dis()\nstandard deviation -> dis(False)\n"+"mathematicaly expected value -> mval()\n"+"probability of deviation from the expectation -> dev()"
            print "calculate probabilities of random value being less than given value -> check()"
            print "build the distribution function of a discrete random variable -> graph()"
        except TypeError , e:
            print e
            print "fill the table with values"
    
    def __check_data(self):
        if type(self.table[self.table.keys()[0]]) not in [type(int()),type(float())] :
            self.error_exeption('Trying to use discrete class function without initialisation.',"not DRV")
    
    @constructor.decor(("mathematicaly expected value is ", "sum(x**2 * p) "))
    def mval(self,f=True):
        self.__check_data()
        if f:
            return sum([k*self.table[k] for k in self.table])
        else:
            return sum([(k**2)*self.table[k] for k in self.table])            
            
    @constructor.decor(("dispersion is ","standard deviation is "))
    def dis(self,f=True):
        m = self.mval._original(self)
        sq = self.mval._original(self,False)
        if f:
            return sq - m**2
        else:
            return round(sqrt(sq - m**2))
    
    @constructor.decor(("probability of deviation from the expectation is ",""))
    def dev(self,f = True):
        m = self.mval._original(self)
        q = self.dis._original(self,False)
        sortable = self.sortedvalues()
        count = self.function_top(sortable)
        return count(m+q)-count(m-q)
    
    @constructor.decor(("probability of random value being less than given value is ",""))
    def check(self,x):
        self.__check_data()
        if type(x) != type(True):
            return self.function_top(self.sortedvalues())(x)
        else:
            return "uncountable, no value was passed"
    
class continuous(constructor):
    def __init__(self,data):
        if type(data) == type(dict()) and type(data.keys()[0]) == type(tuple()): # change [0] condition
            self.table = self.__prep_func(data)
            #self.borders()
        else:
            print "error in initiation, dict expected"
    
    def help(self):
        print "Continuous Random Value"
        print "The function of distribution density -> den()"
        print "The are under the function of distribution density -> area()"
        print "calculate probabilities of random value being less than given value -> check()"
        print "build the distribution function and a function of distribution density of a continuous random variable -> graph()"
        
    def __check_data(self):
        if type(self.table[self.table.keys()[0]]) != type(tuple()):
            self.error_exeption('Trying to use continuous class function without initialisation.',"not CRV")
    
    @constructor.decor(("The function of distribution density is ",""))
    def den(self,f=True):
        self.__check_data()
        res = {}
        for u in self.table:
            res[u] = (self.table[u][0],diff(self.table[u][1],symbols(self.table[u][0])))
        return res
    
    @constructor.decor(("The area under the function of distribution density is ",""))
    def area(self,f=True):
        d_fun = self.den._original(self)
        return sum([integrate(d_fun[u][1],(symbols(d_fun[u][0]),u[0],u[1])) for u in d_fun])
                
    def borders(self):
        return None
    @constructor.decor(("probability of random value being less than given value is",""))
    def check(self,x):
        self.__check_data()
        if type(x) != type(True):
            return self.ice_top(self.sortedvalues())(x)
        else:
            return "uncountable, no value was passed"
        
        
    
    def __prep_func(self,data):
        prep_func = {}
        alp = string.ascii_letters
        for u in data:
            if type(data[u]) == type(str()):
                k = "x"
                for t in data[u]:
                    if t in alp:
                        k = t 
            else:
                k = "const"
            prep_func[(sympify(u[0]),sympify(u[1]))] = (k,sympify(data[u]))
        return prep_func 
            
    
class pack(discrete,continuous):
    def __init__(self,data):
        self.flag = True
        for u in data:
            if type(u) == type(tuple()) :
                self.flag = False
                break
        if self.flag :
            discrete.__init__(self,data)
        else: 
            continuous.__init__(self,data)
            
    def help(self):
        if self.flag:
            discrete.help(self)
        else:
            continuous.help(self)
    
    def check(self,x):
        if self.flag:
            discrete.check(self,x)
        else:
            continuous.check(self,x)
    
    def author(self):
        print "designed and created by\npgrbnv"
    
    
    
    
    
    
