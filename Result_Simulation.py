# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:35:36 2019

@author: 10624
"""

def UnitLength(Link):
    A = Link
    while(len(A) != 20):
        A.append(A[-1])
    return A


LonIter = np.arange(-95, -84.9, 0.1)
LatIter = np.arange(30, 40.1, 0.1)
IntIter = np.arange(0, 10.2, 0.2)

Sys_Water_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sys_Power_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sys_Gas_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sys_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))

Sin_Water_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sin_Power_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sin_Gas_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))
Sin_Sys_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))

Diff_Perform = np.array([[[None]*len(IntIter)]*len(LatIter)]*len(LonIter))

AnalType = 'SingleSum'
Simu_Time = 0

for i in range(len(LonIter)):
    for j in range(len(LatIter)):
        for k in range(len(IntIter)):
            Simu_Time += 1
            print(Simu_Time)
            
            ##System Simulation
            Earth = EarthquakeSys(Shelby_County, LonIter[i], LatIter[j], IntIter[k])
            Earth.DistanceCalculation()
            Earth.NodeFailProbCalculation()
            Earth.MCFailureSimulation()            
            
            while(1):
                Earth.AdjUpdate()
                Earth.FlowUpdate()
                Earth.CascadFail()
                Earth.Performance(AnalType)
                if(Earth.NodeFailIndex[-1] == Earth.NodeFailIndex[-2]):
                    break
            
            Sys_Water_Perform[i][j][k] = np.array(UnitLength(Water.Performance))
            Sys_Power_Perform[i][j][k] = np.array(UnitLength(Water.Performance))
            Sys_Gas_Perform[i][j][k] = np.array(UnitLength(Water.Performance))
            Sys_Perform[i][j][k] = np.array(UnitLength(list(Shelby_County.Performance)))


            ##Single Network Simulation
            for Network in Shelby_County.Networks:
                Network.FlowAdj = []
                Network.FlowAdj.append(np.zeros([Network.NodeNum, Network.NodeNum]))
                m = Network.WholeNodeSeries[Network.SupplySeries[0]]
                n = Network.WholeNodeSeries[Network.DemandSeries[-1]]
                Network.FlowAdj[0] = Shelby_County.FlowAdj[0][m:(n + 1), m:(n + 1)]
            
            for Network in Shelby_County.Networks:
                exec('Earth{} = EarthquakeNet(Network, Earth)'.format(Network.Name))
                exec('Earth{}.VariableIni()'.format(Network.Name))
                exec('Earth{}.IniFailCopy()'.format(Network.Name))
                while(1):
                    exec('Earth{}.AdjUpdate()'.format(Network.Name))
                    exec('Earth{}.FlowUpdate()'.format(Network.Name))
                    exec('Earth{}.CascadFail()'.format(Network.Name))
                    exec('Earth{}.Performance("{}")'.format(Network.Name, AnalType))
                    exec('Temp1 = Earth{}.Target.NodeFailIndex[-1]'.format(Network.Name))
                    exec('Temp2 = Earth{}.Target.NodeFailIndex[-2]'.format(Network.Name))
                    if(Temp1 == Temp2):
                        break
            
            Sin_Water_Perform[i][j][k] = np.array(UnitLength(Water.SinglePerform))
            Sin_Power_Perform[i][j][k] = np.array(UnitLength(Power.SinglePerform))
            Sin_Gas_Perform[i][j][k] = np.array(UnitLength(Gas.SinglePerform))
            Sin_Sys_Perform[i][j][k] = (Sin_Water_Perform[i][j][k] + Sin_Power_Perform[i][j][k] + Sin_Gas_Perform[i][j][k])/3
            
            Diff_Perform[i][j][k] = (Sin_Sys_Perform[i][j][k] - Sys_Perform[i][j][k])/Sys_Perform[i][j][k]