from netCDF4 import Dataset
import numpy as np
import math

#-------------------------------------------------------------

def wind_velocity(x, y, Lx, Ly):
    
    dc = 1000.0
    #mx = dc*256.0
    #my = dc*256.0
    mx = 256000.0 #+ 51200.0*time/86400.0
    my = 256000.0 #+ 51200.0*time/86400.0
    
    r = math.sqrt((mx-x)**2 + (my-y)**2)
    alpha = 72.0
    
    svmax = 3.0e-4*math.exp(-1e-5*r)
    #svmax = 0.3*math.exp(-r/16e5)
    u = -svmax*(math.cos(alpha)*(x-mx) + math.sin(alpha)*(y-my))#/Lx
    v = -svmax*(-math.sin(alpha)*(x-mx) + math.cos(alpha)*(y-my))#/Ly
    
    #u = 5.0  - 3.0 * math.sin((2.0 * math.pi * x) / Lx) * math.sin((math.pi * y) / Ly)
    #v = 5.0  - 3.0 * math.sin((2.0 * math.pi * y) / Ly) * math.sin((math.pi * x) / Lx)

    return u, v

#-------------------------------------------------------------

def ocean_currents(x, y, Lx, Ly):

    u =  0.01 * ((2.0 * y - Ly) / Ly)
    v = -0.01 * ((2.0 * x - Ly) / Ly)

    return u, v

#-------------------------------------------------------------

def ice_concentration(x, y, Lx, Ly):

    conc = 1.0
    #conc = max(min(x / Lx, 1.0),0.0)
    
    return conc

#-------------------------------------------------------------

def ice_thickness(x,y):

    dc = 16000.0
    H = 0.3 + 0.005*(math.sin((60.0/1.0e6)*x) + math.sin((30.0/1.0e6)*y))

    return H
#-------------------------------------------------------------

def create_ic(gridfile, icfile):

    # load grid file
    grid = Dataset(gridfile, "r")

    xCell = grid.variables["xCell"][:]
    yCell = grid.variables["yCell"][:]

    xVertex = grid.variables["xVertex"][:]
    yVertex = grid.variables["yVertex"][:]

    nCells = len(grid.dimensions["nCells"])
    nVertices = len(grid.dimensions["nVertices"])
    nEdges = len(grid.dimensions["nEdges"])

    # calculate output variables
    uAirVelocity = np.empty(nCells)
    vAirVelocity = np.empty(nCells)

    uOceanVelocity = np.empty(nCells)
    vOceanVelocity = np.empty(nCells)

    iceConcentration = np.empty((nCells,1,1))
    iceVolume        = np.empty((nCells,1,1))

    fVertex = np.empty((nVertices,1,1))
    fEdge = np.empty((nEdges,1,1))

    xmin = np.amin(xVertex)
    xmax = np.amax(xVertex)
    Lx = xmax - xmin

    print(xmin, xmax, Lx)

    ymin = np.amin(yVertex)
    ymax = np.amax(yVertex)
    Ly = ymax - ymin

    print(ymin, ymax, Ly)
    L = 512000.0

    for iCell in range(0,nCells):

        x = xCell[iCell] - xmin
        y = yCell[iCell] - ymin

        uAirVelocity[iCell],   vAirVelocity[iCell]   = wind_velocity (x, y, Lx, Ly)
        uOceanVelocity[iCell], vOceanVelocity[iCell] = ocean_currents(x, y, Lx, Ly)

        iceConcentration[iCell,0,0] = ice_concentration(x, y, Lx, Ly)

        iceVolume[iCell,0,0] = ice_thickness(x,y) * iceConcentration[iCell,0,0]
        #iceVolume[iCell,0,0] = 2.0 * iceConcentration[iCell,0,0]

    for iVertex in range(0,nVertices):
        fVertex[iVertex] = 1.46e-4
    for iEdge in range(0,nEdges):
        fEdge[iEdge] = 1.46e-4

    # create output file
    ic = Dataset(icfile, "w", format="NETCDF3_64BIT")
    # format: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA

    ic.createDimension("nCells", nCells)
    ic.createDimension("nVertices", nVertices)
    ic.createDimension("nEdges", nEdges)
    ic.createDimension("nCategories", 1)
    ic.createDimension("ONE", 1)
    ic.createDimension("Time", None)

    uAirVelocityVar = ic.createVariable("uAirVelocity", 'd', ('nCells'))
    vAirVelocityVar = ic.createVariable("vAirVelocity", 'd', ('nCells'))
    uAirVelocityVar[:] = uAirVelocity[:]
    vAirVelocityVar[:] = vAirVelocity[:]

    uOceanVelocityVar = ic.createVariable("uOceanVelocity", 'd', ('nCells'))
    vOceanVelocityVar = ic.createVariable("vOceanVelocity", 'd', ('nCells'))
    uOceanVelocityVar[:] = uOceanVelocity[:]
    vOceanVelocityVar[:] = vOceanVelocity[:]

    iceConcentrationVar = ic.createVariable("iceAreaCategory", 'd', ('nCells', 'nCategories', 'ONE'))
    iceConcentrationVar[:,:,:] = iceConcentration[:,:,:]

    iceVolumeVar = ic.createVariable("iceVolumeCategory", 'd', ('nCells', 'nCategories', 'ONE'))
    iceVolumeVar[:,:,:] = iceVolume[:,:,:]

    fVertexVar = ic.createVariable("fVertex", 'd', ('nVertices'))
    fVertexVar[:] = fVertex[:]
    fEdgeVar = ic.createVariable("fEdge", 'd', ('nEdges'))
    fEdgeVar[:] = fEdge[:]

    ic.close()
    grid.close()

#-------------------------------------------------------------

def create_ics():

    gridTypes = ["hex"]
    #gridTypes = ["quad"]
    #grids = {"hex": ["0082x0094",
    #                 "0164x0188",
    #                 "0328x0376",
    #                 "0656x0752"],
    #         "quad":["0080x0080",
    #                 "0160x0160",
    #                 "0320x0320",
    #                 "0640x0640"]}

    #grids = {"hex": ["0514x0526"],
    #         "quad":["0512x0512"]}
    grids = {"hex": ["2562"]}

    for gridType in gridTypes:
        for grid in grids[gridType]:

            gridfile = "grid_%s_%s.nc" %(gridType,grid)
            icfile = "ic_%s_%s.nc" %(gridType,grid)

            create_ic(gridfile, icfile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics()
