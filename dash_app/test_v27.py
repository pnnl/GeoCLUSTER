#precalculate the thermal response with a line and cylindrical heat source. Precalculating allows to speed up the SBT algorithm.
#precalculate finite pipe correction
fpcminarg = min(Deltaz)**2 / (4 * alpha_m * times[-1])
fpcmaxarg = max(Deltaz)**2 / (4 * alpha_m * (min(times[1:] - times[:-1])))
Amin1vector = np.logspace(np.log10(fpcminarg) - 0.1, np.log10(fpcmaxarg) + 0.1, NoArgumentsFinitePipeCorrection)
finitecorrectiony = np.zeros(NoArgumentsFinitePipeCorrection)

for i, Amin1 in enumerate(Amin1vector):
    Amax1 = (16)**2
    if Amin1 > Amax1:
        Amax1 = 10 * Amin1
    Adomain1 = np.logspace(np.log10(Amin1), np.log10(Amax1), NoDiscrFinitePipeCorrection)
    finitecorrectiony[i] = np.trapz(-1 / (Adomain1 * 4 * np.pi * k_m) * erfc(1/2 * np.power(Adomain1, 1/2)), Adomain1)

#precalculate besselintegration for infinite cylinder
if clg_configuration == 1: #co-axial geometry
    besselminarg = alpha_m * (min(times[1:] - times[:-1])) / radius**2
    besselmaxarg = alpha_m * timeforlinesource / radius**2
elif clg_configuration == 2: #U-loop geometry
    besselminarg = alpha_m * (min(times[1:] - times[:-1])) / max(radiusvector)**2
    besselmaxarg = alpha_m * timeforlinesource / min(radiusvector)**2
deltazbessel = np.logspace(-10, 8, NoDiscrInfCylIntegration)
argumentbesselvec = np.logspace(np.log10(besselminarg) - 0.5, np.log10(besselmaxarg) + 0.5, NoArgumentsInfCylIntegration)
besselcylinderresult = np.zeros(NoArgumentsInfCylIntegration)

for i, argumentbessel in enumerate(argumentbesselvec):
    besselcylinderresult[i] = 2 / (k_m * np.pi**3) * np.trapz((1 - np.exp(-deltazbessel**2 * argumentbessel)) / (deltazbessel**3 * (jv(1, deltazbessel)**2 + yv(1, deltazbessel)**2)), deltazbessel)
print("SBT distributions calculated successfully")

N = len(Deltaz)  # Number of elements
elementcenters = 0.5 * np.column_stack((x[1:], y[1:], z[1:])) + 0.5 * np.column_stack((x[:-1], y[:-1], z[:-1]))  # Matrix that stores the mid point coordinates of each element
if clg_configuration == 2: #U-loop geometry
    interconnections = interconnections - 1
    elementcenters = np.delete(elementcenters, interconnections.reshape(-1,1), axis=0)  # Remove duplicate coordinates

SMatrix = np.zeros((N, N))  # Initializes the spacing matrix, which holds the distance between center points of each element [m]
SoverL = np.zeros((N, N))  # Initializes the ratio of spacing to element length matrix
for i in range(N):
    SMatrix[i, :] = np.sqrt((elementcenters[i, 0] - elementcenters[:, 0])**2 + (elementcenters[i, 1] - elementcenters[:, 1])**2 + (elementcenters[i, 2] - elementcenters[:, 2])**2)
    SoverL[i, :] = SMatrix[i, :] / Deltaz[i] #Calculates the ratio of spacing between two elements and element length

