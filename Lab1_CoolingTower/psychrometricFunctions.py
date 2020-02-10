# Fig 1c
import psychrolib as psy
import numpy as np
# Set the unit system, for example to SI (can be either psychrolib.SI or psychrolib.IP)
psy.SetUnitSystem(psy.SI)

amb_pressure = 86846.203 #pa
Tdb = np.array([[20.5, 19.8, 18.1, 20.2, 21.5], [20.8, 19.9, 17.4, 19.5, 21.4], [21.2, 20.2, 19.4, 21.2, 21.9]])
Twb = np.array([[10.0, 13.9, 17.2, 19.8, 19.0], [10.3, 13.1, 16.2, 18.6, 18.7], [10.3, 15.6, 18.8, 21.1, 19.6]])
spec_vol_air = np.zeros([Tdb.shape[0], Tdb.shape[1]])
spec_hum = np.zeros([Tdb.shape[0], Tdb.shape[1]])

for i in range(Tdb.shape[0]):
    for j in range(Tdb.shape[1]):
        spec_hum[i, j] = psy.GetSpecificHumFromHumRatio(psy.GetHumRatioFromTWetBulb(Tdb[i,j], Twb[i,j], amb_pressure))
        # spec_vol_air[i,j] = psy.GetDryAirVolume(Tdb[i,j], amb_pressure)
        spec_vol_air[i,j] = psy.GetMoistAirVolume(Tdb[i, j], psy.GetHumRatioFromTWetBulb(Tdb[i,j], Twb[i,j], amb_pressure), amb_pressure)

print(spec_vol_air)

# fig 1e
print(spec_hum)


# fig 1f
for i in range(3):
    enthalpy_in = psy.GetDryAirEnthalpy(Tdb[i,0])
    enthalpy_out = psy.GetDryAirEnthalpy(Tdb[i,-1])
    print('mass idx, in, out: ',i, '\t', enthalpy_in, '\t', enthalpy_out)


