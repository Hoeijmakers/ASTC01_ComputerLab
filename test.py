import TransitSimulator

sys1 = System('WASP-121b', 1.183, 1.865, 3.34, 0.0558, 1.45, 1.5, 0,  np.pi / 2., 0)

sys1.TransitParams()
sys1.Transit()
