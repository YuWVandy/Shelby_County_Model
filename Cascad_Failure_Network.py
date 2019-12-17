# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 16:56:24 2019

@author: wany105
"""


class EarthquakeNet:
    def __init__(self, Target, FailureSystem):##Type 1 - System, Type 2 - Single Network
        self.Target = Target
        self.FailSys = FailureSystem
    
    def VariableIni(self):
        self.Target.SingleTimeAdj = []
        self.Target.SingleTimeAdj.append(copy.copy(Network.Adj))
        self.Target.NodeFailIndex = []
        self.Target.SingleSatisDe = []
        self.Target.SinglePerform = []

    def IniFailCopy(self):
        for node in self.Target.NodeSeries:
            if(self.Target.WholeNodeSeries[node] in self.FailSys.NodeFailIndex[0]):
                self.Target.NodeFailIndex.append(node)
        
    def AdjUpdate(self):
        Adj = copy.copy(self.Target.TimeAdj[-1])
        
        Adj[self.Target.NodeFailIndex[-1], :] = 0
        Adj[:, self.Target.NodeFailIndex[-1]] = 0

        self.Target.TimeAdj.append(Adj)
        
    def FlowUpdate(self):
        Flow = self.Target.TimeAdj[-1]*self.Target.FlowAdj[-1]
        ##Demand -> Supply
        for node in self.Target.SupplySeries:
            FlowInNode = np.sum(Flow[:, node])
            Ratio = Flow[node, :]/np.sum(Flow[node, :])
            Ratio[np.isnan(Ratio)] = 0
            if(np.sum(Flow[node, :]) != 0):
                Flow[node, :] = FlowInNode*Ratio
        ##Supply -> Demand
        self.Target.SingleSatisDe.append([])
        for node in self.Target.DemandSeries:
            FlowInNode = np.sum(Flow[self.Target.SupplySeries, node]) + np.sum(Flow[self.Target.DemandSeries, node])
            if(FlowInNode >= self.Target.DemandValue[node - self.Target.DemandSeries[0]]):
                FlowInNode = FlowInNode - self.Target.DemandValue[node - self.Target.DemandSeries[0]]
                self.Target.SingleSatisDe[-1].append(self.Target.DemandValue[node - self.Target.DemandSeries[0]])
            else:
                self.Target.SingleSatisDe[-1].append(0.25*FlowInNode)
                FlowInNode *= 0.75
                
            Ratio = Flow[node, :]/np.sum(Flow[node, :])
            Ratio[np.isnan(Ratio)] = 0
            if(np.sum(Flow[node, :]) != 0):
                Flow[node, :] = FlowInNode*Ratio        
                
        self.Target.FlowAdj.append(Flow)
        
    def CascadFail(self):
        print(self.Target.NodeFailIndex)
        self.Target.NodeFailIndex.append(copy.copy(self.Target.NodeFailIndex[-1]))
        for i in range(Network.NodeNum):
            FlowInNode = np.sum(self.Target.FlowAdj[-1][:, i])
            FlowInNode0 = np.sum(self.Target.FlowAdj[0][:, i])
            
            if(FlowInNode < LowBound*FlowInNode0 and (i not in self.Target.NodeFailIndex[-1])):
                self.Target.NodeFailIndex[-1].append(i)
            if(FlowInNode > UpBound*FlowInNode0 and (i not in self.Target.NodeFailIndex[-1])):
                self.Target.NodeFailIndex[-1].append(i)

    def Performance(self, Type):
        self.Target.Performance.append([])
        if(Type == "SingleSum"):
            Temp = np.array(self.Target.SingleSatisDe[-1])/self.Target.DemandValue
            Temp[np.isnan(Temp)] = 0
            for i in range(self.Target.DemandNum):
                if(Temp[i] > 1):
                    Temp = 1
            self.Target.Performance[-1] = min(1, np.average(Temp))
        
        if(Type == "WholeSum"):
            self.Target.Performance[-1] = min(1, np.sum(np.array(self.Target.SingleSatisDe[-1]))/np.sum(np.array(self.Target.SatisfyDemand[0])))
        
for Network in Shelby_County.Networks:
    exec('Earth{} = EarthquakeNet(Network, Earth)'.format(Network.Name))
    exec('Earth{}.VariableIni()'.format(Network.Name))
    exec('Earth{}.IniFailCopy()'.format(Network.Name))
    while(1):
        exec('Earth{}.AdjUpdate()'.format(Network.Name))
        exec('Earth{}.FlowUpdate()'.format(Network.Name))
        exec('Earth{}.CascadFail()'.format(Network.Name))
        exec('Earth{}.Performance("{}")'.format(Network.Name, Type))
        exec('Temp1 = Earth{}.Target.NodeFailIndex[-1]'.format(Network.Name))
        exec('Temp2 = Earth{}.Target.NodeFailIndex[-2]'.format(Network.Name))
        if(Temp1 == Temp2):
            break