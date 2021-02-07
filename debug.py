import importlib.util
spec = importlib.util.spec_from_file_location("hurdat", "../ORB/python/convective-structure/hurdat.py")
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)

hdata = h.Hurdat()

a = hdata.genesis_to_lysis_filter(50)

a.identify_events(threshold = 25)