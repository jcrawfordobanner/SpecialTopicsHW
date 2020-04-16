import picos as pic
from picos import RealVariable
from copy import deepcopy
from heapq import *
import heapq as hq
import numpy as np
import itertools
import math
counter = itertools.count()

class BBTreeNode():
    def __init__(self, vars = [], constraints = [], objective='', prob=None):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.prob = prob

    def __deepcopy__(self, memo):
        '''
        Deepcopies the picos problem
        This overrides the system's deepcopy method bc it doesn't work on classes by itself
        '''
        newprob = pic.Problem.clone(self.prob)
        return BBTreeNode(self.vars, newprob.constraints, self.objective, newprob)

    def buildProblem(self):
        '''
        Bulids the initial Picos problem
        '''
        prob=pic.Problem()

        prob.add_list_of_constraints(self.constraints)

        prob.set_objective('max', self.objective)
        self.prob = prob
        return self.prob

    def is_integral(self):
        '''
        Checks if all variables (excluding the one we're maxing) are integers
        '''
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
                return False
        return True

    def branch_floor(self, branch_var):
        '''
        Makes a child where xi <= floor(xi)
        '''
        n1 = deepcopy(self)
        n1.prob.add_constraint( branch_var <= math.floor(branch_var.value) ) # add in the new binary constraint

        return n1

    def branch_ceil(self, branch_var):
        '''
        Makes a child where xi >= ceiling(xi)
        '''
        n2 = deepcopy(self)
        n2.prob.add_constraint( branch_var >= math.ceil(branch_var.value) ) # add in the new binary constraint
        return n2

    def integral(self,v):
        '''
        Checks if all variables (excluding the one we're maxing) are integers
        '''
        if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
            return False
        return True

    def bbsolve(self):
        '''
        Use the branch and bound method to solve an integer program
        This function should return:
            return bestres, bestnode_vars

        where bestres = value of the maximized objective function
              bestnode_vars = the list of variables that create bestres
        '''

        # these lines build up the initial problem and adds it to a heap
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        heap = [(res, next(counter), root)]
        bestres = -1e20 # a small arbitrary initial best objective value
        bestnode_vars = root.vars # initialize bestnode_vars to the root vars
        #TODO: fill this part in
        best_prob=root.prob
        best_sol=res
        feasible=True
        int1=False
        int2=False
        current=root
        while(feasible):
            if(current.is_integral()):
                best_sol.apply(toProblem=best_prob)
                return bestres, bestnode_vars
            else:
                feasible1=False
                feasible2=False
                for i in self.vars[:-1]:
                    if(not current.integral(i)):
                        feasible1=False
                        feasible2=False
                        branch1=current.branch_ceil(i)
                        try:
                            sol1=branch1.prob.solve(solver='cvxopt')
                        except:
                            sol1=False
                        if(sol1!=False and len(sol1.primals)!=0):
                            feasible1=True
                            int1=branch1.is_integral()
                        res.apply(toProblem=current.prob)
                        branch2=current.branch_floor(i)
                        try:
                            sol2=branch2.prob.solve(solver='cvxopt')
                        except:
                            sol2=False
                        if(sol2!=False and len(sol2.primals)!=0):
                            feasible2=True
                            int2=branch2.is_integral()
                        if(feasible1 and feasible2):
                            sol1.apply(toProblem=branch1.prob)
                            temp1=float(branch1.vars[-1])
                            sol2.apply(toProblem=branch2.prob)
                            temp2=float(branch2.vars[-1])
                            if(int1 and int2):
                                if(temp1>=temp2):
                                    sol1.apply(toProblem=branch1.prob)
                                    if(float(branch1.vars[-1])>bestres):
                                        bestres=float(branch1.vars[-1])
                                        bestnode_vars=branch1.vars
                                        best_prob=branch1.prob
                                        best_sol=sol1
                                    best_sol.apply(toProblem=best_prob)
                                    return bestres, bestnode_vars
                                elif(temp1<temp2):
                                    sol2.apply(toProblem=branch2.prob)
                                    if(float(branch2.vars[-1])>bestres):
                                        bestres=float(branch2.vars[-1])
                                        bestnode_vars=branch2.vars
                                        best_prob=branch2.prob
                                        best_sol=sol2
                                    best_sol.apply(toProblem=best_prob)
                                    return bestres, bestnode_vars
                            elif(int1):
                                if(temp1>=temp2):
                                    sol1.apply(toProblem=branch1.prob)
                                    if(float(branch1.vars[-1])>bestres):
                                        bestres=float(branch1.vars[-1])
                                        bestnode_vars=branch1.vars
                                        best_prob=branch1.prob
                                        best_sol=sol1
                                    best_sol.apply(toProblem=best_prob)
                                    return bestres, bestnode_vars
                                elif(temp1<temp2):
                                    sol1.apply(toProblem=branch1.prob)
                                    if(float(branch2.vars[-1])>bestres):
                                        bestres=float(branch1.vars[-1])
                                        bestnode_vars=branch1.vars
                                        best_prob=branch1.prob
                                        best_sol=sol1
                                    sol2.apply(toProblem=branch2.prob)
                                    current=branch2
                                    res=sol2
                                    break
                            elif(int2):
                                if(temp1>=temp2):
                                    sol2.apply(toProblem=branch2.prob)
                                    if(float(branch1.vars[-1])>bestres):
                                        bestres=float(branch2.vars[-1])
                                        bestnode_vars=branch2.vars
                                        best_prob=branch2.prob
                                        best_sol=sol2
                                    sol1.apply(toProblem=branch1.prob)
                                    current=branch1
                                    res=sol1
                                elif(temp1<temp2):
                                    sol2.apply(toProblem=branch2.prob)
                                    if(float(branch2.vars[-1])>bestres):
                                        bestres=float(branch2.vars[-1])
                                        bestnode_vars=branch2.vars
                                        best_prob=branch2.prob
                                        best_sol=sol2
                                    best_sol.apply(toProblem=best_prob)
                                    return bestres, bestnode_vars
                            else:
                                if(temp1>=temp2):
                                    sol1.apply(toProblem=branch1.prob)
                                    if(float(branch1.vars[-1])>bestres):
                                        bestres=float(branch1.vars[-1])
                                        bestnode_vars=branch1.vars
                                        best_prob=branch1.prob
                                        best_sol=sol1
                                    current=branch1
                                    res=sol1
                                elif(temp1<temp2):
                                    sol2.apply(toProblem=branch2.prob)
                                    if(float(branch2.vars[-1])>bestres):
                                        bestres=float(branch2.vars[-1])
                                        bestnode_vars=branch2.vars
                                        best_prob=branch2.prob
                                        best_sol=sol2
                                    current=branch2
                                    res=sol2
                        elif(feasible1):
                            sol1.apply(toProblem=branch1.prob)
                            if(float(branch1.vars[-1])>bestres and branch1.is_integral()):
                                bestres=float(branch1.vars[-1])
                                bestnode_vars=branch1.vars
                                best_prob=branch1.prob
                                best_sol=sol1
                            current=branch1
                            res=sol1
                            break
                        elif(feasible2):
                            sol2.apply(toProblem=branch2.prob)
                            if(float(branch2.vars[-1])>bestres and branch2.is_integral()):
                                bestres=float(branch2.vars[-1])
                                bestnode_vars=branch2.vars
                                best_prob=branch2.prob
                                best_sol=sol2
                            current=branch2
                            res=sol2
                            break
                        else:
                            feasible=False
        best_sol.apply(toProblem=best_prob)
        return bestres, bestnode_vars
