
#rootcint -f BeamDictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p Beam.h BeamLinkDef.h
#g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c BeamDictionary.cpp -o BeamDictionary.o
#g++ -g -fPIC -shared -o/home/phong/lib/libBeam.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include Beam.cpp BeamDictionary.cpp

#rootcint -f CloverDictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p Clover.h CloverLinkDef.h
#g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c CloverDictionary.cpp -o CloverDictionary.o
#g++ -g -fPIC -shared -o/home/phong/lib/libClover.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include Clover.cpp CloverDictionary.cpp

#rootcint -f BELENDictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p BELEN.h BELENLinkDef.h
#g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c BELENDictionary.cpp -o BELENDictionary.o
#g++ -g -fPIC -shared -o/home/phong/lib/libBELEN.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include BELEN.cpp BELENDictionary.cpp

#rootcint -f AIDADictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p AIDA.h AIDALinkDef.h
#g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c AIDADictionary.cpp -o AIDADictionary.o
#g++ -g -fPIC -shared -o/home/phong/lib/libAIDA.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include AIDA.cpp AIDADictionary.cpp

rootcint -f TreeDataDictionary.cpp -c -Wall -fPIC -O3 `root-config --cflags` -p TreeData.h TreeDataLinkDef.h
#g++ -Wall -Wno-long-long -g -O3 `root-config --cflags` -fPIC -c AIDADictionary.cpp -o AIDADictionary.o
#g++ -g -fPIC -shared -o/home/phong/lib/libAIDA.so `root-config --ldflags` -Wall -fPIC -O3 -I/home/phong/root/include AIDA.cpp AIDADictionary.cpp

