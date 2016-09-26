rootcint -f BeamDictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p Beam.h BeamLinkDef.h
g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c BeamDictionary.cpp -o BeamDictionary.o
g++ -g -fPIC -shared -o/home/phong/lib/libBeam.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include Beam.cpp BeamDictionary.cpp

rootcint -f AIDADictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p AIDA.h AIDALinkDef.h
g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c AIDADictionary.cpp -o AIDADictionary.o
g++ -g -fPIC -shared -o/home/phong/lib/libAIDA.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include AIDA.cpp AIDADictionary.cpp

rootcint -f EURICADictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p EURICA.h EURICALinkDef.h
g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c EURICADictionary.cpp -o EURICADictionary.o
g++ -g -fPIC -shared -o/home/phong/lib/libEURICA.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include EURICA.cpp EURICADictionary.cpp




