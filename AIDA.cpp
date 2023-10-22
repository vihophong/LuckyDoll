#include "AIDA.h"

bool AIDA::BetaGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;
                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());
                cluster->SetRankingFlag(1);

                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    //! Deallocating memory
                    delete cluster;
                }
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;

}


bool AIDA::BetaGetPosAllNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< multimap<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    multimap<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;
                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());


                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());



                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        ConstructEdiffRankingFlags();
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}





bool AIDA::BetaGetPosYonly()
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short y_prev = -11;
        short nStripInClusterY = 0;
        double E_Y = 0;
        double E_Y_ch = 0;
        nClusterY[z] = 0;//ok...
        //! loop through the  y axis
        for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
            short y = yhitmaps_it->first;
            double yen = yhitmaps_it->second->GetEnergy();
            //!out of cluster
            if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                if (nStripInClusterY<=maxStrInCluster){
                    posY = E_Y_ch/E_Y;
                    //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                    AIDACluster* cluster = new AIDACluster;
                    cluster->SetHitPosition(0,posY,z);
                    cluster->SetXEnergy(0);
                    cluster->SetYEnergy(E_Y);
                    cluster->SetXMult(0);
                    cluster->SetYMult(nStripInClusterY);
                    //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                    //! Timming
                    short posYint = (short) round(posY);
                    AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;
                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());

                    //! finally insert map
                    cluster_map[z].insert(make_pair(-E_Y,cluster));
                    nClusterY[z]++;
                }
                //! reset hits and energy in Y
                E_Y = 0;
                E_Y_ch = 0;
                nStripInClusterY = 0;
            }

            //!still in a same cluster
            if (yen>0){
                //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                E_Y += yen;
                E_Y_ch += (double) y * yen;
                nStripInClusterY++;
                //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
            }

            //! handle last cluster!
            if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                if (nStripInClusterY<=maxStrInCluster){
                    posY = E_Y_ch/E_Y;
                    //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                    AIDACluster* cluster = new AIDACluster;
                    cluster->SetHitPosition(0,posY,z);
                    cluster->SetXEnergy(0);
                    cluster->SetYEnergy(E_Y);
                    cluster->SetXMult(0);
                    cluster->SetYMult(nStripInClusterY);

                    //! Timming
                    short posYint = (short) round(posY);
                    AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;
                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());


                    //! if there is available fast time stamp
                    //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                    //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                    //! finally insert map
                    cluster_map[z].insert(make_pair(-E_Y,cluster));
                    nClusterY[z]++;
                }
                //! reset hits and energy in Y
                E_Y = 0;
                E_Y_ch = 0;
                nStripInClusterY = 0;
            }
            y_prev = y;
        }//y
        fnyclustersz[z]=nClusterY[z];
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        if (nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            Int_t icluster=0;
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());
                cluster->SetRankingFlag(icluster);
                this->AddCluster(cluster);
                maxZ = z;
                icluster++;
            }//loop all cluster map

        }//if condition
    }//all z

    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}




bool AIDA::BetaGetPosAllNewMax(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = -1;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch;//max
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = -1;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            if (y>E_Y_ch) E_Y_ch = (double) y;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = -1;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                if (x>E_X_ch) E_X_ch = (double) x;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch;//maximum
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = -1;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            if (y>E_Y_ch) E_Y_ch=(double) y;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());


                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = -1;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        fnyclustersz[z]=nClusterY[z];
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}




bool AIDA::IonGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        fnyclustersz[z]=nClusterY[z];
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;
                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());
                cluster->SetRankingFlag(1);

                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    //! Deallocating memory
                    delete cluster;
                }
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;

}


bool AIDA::IonGetPosAllNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
            }

            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch/E_Y;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());


                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z




    if(this->GetNClusters()>0) {
        ConstructEdiffRankingFlags();
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}


bool AIDA::IonGetPosAllNewMax(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
            //if (z==0&&hit->GetEnergy()>0) cout<<xy<<"b-"<<hit->GetEnergy()<<endl;
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    map<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = -1;
        //! loop through the x axis

        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            //if (this->GetHit(0)->GetTimestamp()==23112644336) cout<<z<<"eee-"<<x<<"-"<<xen<<endl;
            //if (z==0) cout<<x<<"cc-"<<xen<<endl;

            //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"*"<<endl;

            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch;//max
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = -1;
                    nClusterY[z] = 0;//ok...
                    //! loop through the  y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){

                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                //cout<<"x"<<posX <<"y"<<posY<<"e" <<(E_Y/E_X-1)*(E_Y/E_X-1)<<endl;

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            if (y>E_Y_ch) E_Y_ch = (double) y;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        y_prev = y;
                    }//y
                     nClusterX[z]++;
                }//copy to the last cluster handling in X

                //****************************************//

                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = -1;
                nStripInClusterX = 0;
            }
            //!still in a same cluster
            if ( xen>0){
                //if (thereis) cout<<"z"<<z<<"x"<<x<<"en"<<xen<<"xprev"<<x_prev<<endl;
                E_X += xen;
                if (x>E_X_ch) E_X_ch = (double) x;
                nStripInClusterX++;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<x<<"dd-"<<xen<<"-"<<nStripInClusterX<<endl;
            }
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                //cout<<"ggg"<<x<<"-"<<E_X_ch/E_X<<endl;
                //if (nClusterX[z]>1&&nStripInClusterX>1) cout<<"\n\n\neee"<<x<<"-"<<E_X_ch/E_X<<endl;

                //****************************************//

                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch;//maximum
                    //if (nStripInClusterX>1) cout<<posX<<"a-"<<E_X_ch<<"-"<<E_X<<endl;

                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = -1;
                    nClusterY[z] = 0;//ok...
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }

                        //!still in a same cluster
                        if (yen>0){
                            //if (thereis) cout<<"z"<<z<<"y"<<y<<"en"<<yen<<"xprev"<<y_prev<<endl;
                            E_Y += yen;
                            if (y>E_Y_ch) E_Y_ch=(double) y;
                            nStripInClusterY++;
                            //if (nStripInClusterY>1) cout<<y<<"dd-"<<yen<<"-"<<nStripInClusterY<<endl;
                        }

                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            //if (thereis)cout<<"eee"<<y<<"-"<<E_Y_ch/E_Y<<endl;
                            if (nStripInClusterY<=maxStrInCluster){
                                posY = E_Y_ch;
                                //if (nStripInClusterY>1&&nClusterY[z]>1) cout<<posY<<"b-"<<E_Y_ch<<"-"<<E_Y<<endl;
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);

                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                /*
                                //! Take the ealiest time stamp
                                if (hitAtX->GetTimestamp() > hitAtY->GetTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                //! if there is available fast time stamp
                                //if (hitAtX->GetFastTimestamp()==0 || hitAtY->GetFastTimestamp()==0)
                                //    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp() + hitAtY->GetFastTimestamp());


                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));
                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = -1;
                            nStripInClusterY = 0;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }//copy to the last cluster handling in X
                //****************************************//
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = -1;
                nStripInClusterX = 0;

            }
            x_prev = x;
        }//x
        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z

    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}

bool AIDA::ConstructEdiffRankingFlags()
{

    //! energy difference EY-EX
    //cluster map sort by (EX-EY)(EX-EY)
    vector< map<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    map<Double_t,AIDACluster*>::iterator cluster_map_it;
    //cluster map sort by EX+EY
    vector< map<Double_t,AIDACluster*> > cluster_map2(NumDSSD);
    map<Double_t,AIDACluster*>::iterator cluster_map2_it;

    //cluster map x y sort by strip number
    std::map < Double_t ,Double_t> fxmap[NumDSSD];
    std::map < Double_t ,Double_t> fymap[NumDSSD];
    std::map < Double_t ,Double_t> ::iterator fxmap_it[NumDSSD];
    std::map < Double_t ,Double_t> ::iterator fymap_it[NumDSSD];


    for (int i=0;i<fnclusters;i++){
        AIDACluster* cluster=this->GetCluster(i);
        Double_t xpos=cluster->GetHitPositionX();
        Double_t ypos=cluster->GetHitPositionY();
        Int_t z=(Int_t)cluster->GetHitPositionZ();
        Double_t ex=cluster->GetXEnergy();
        Double_t ey=cluster->GetYEnergy();
        cluster_map[z].insert(make_pair((ey-ex)*(ey-ex),cluster));
        cluster_map2[z].insert(make_pair(-ey-ex,cluster));

        fxmap[z].insert(std::make_pair(xpos,-1.));
        fymap[z].insert(std::make_pair(ypos,-1.));
    }

    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
            AIDACluster* cluster = cluster_map_it->second;
            if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))
            {
                cluster->SetEDiffRankingFlag(1);
                xindex.push_back(cluster->GetHitPositionX());
                yindex.push_back(cluster->GetHitPositionY());
            }
        }//loop all cluster map

        unsigned short sumexyrank=0;
        for(cluster_map2_it = cluster_map2[z].begin(); cluster_map2_it != cluster_map2[z].end(); cluster_map2_it++){
            AIDACluster* cluster = cluster_map2_it->second;
            cluster->SetSumEXYRank(sumexyrank);
            sumexyrank++;
        }//loop all cluster map

        double xprev=-1;
        double yprev=-1;
        fmindx[z]=9999;
        fmindy[z]=9999;
        for (fxmap_it[z]=fxmap[z].begin();fxmap_it[z]!=fxmap[z].end();fxmap_it[z]++){
            if (xprev>0&&((fxmap_it[z]->first-xprev)<fmindx[z])) fmindx[z]=fxmap_it[z]->first-xprev;
            xprev=fxmap_it[z]->first;
        }
        for (fymap_it[z]=fymap[z].begin();fymap_it[z]!=fymap[z].end();fymap_it[z]++){
            if (yprev>0&&((fymap_it[z]->first-yprev)<fmindy[z])) fmindy[z]=fymap_it[z]->first-yprev;
            yprev=fymap_it[z]->first;
        }
    }//all z
    return true;
}


bool AIDA::BetaGetPosAllNew2(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< multimap<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    multimap<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){

        //! start with X clusters
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        unsigned long long tsminx=0;
        unsigned long long tsmaxx=0;
        short minPosX=0;
        short maxPosX=0;
        double MAX_E_X = -9999999;
        short X_MAX_E_X = -1;

        //! loop through the x axis
        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            unsigned long long xts = xhitmaps_it->second->GetTimestamp();
            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    maxPosX = x_prev;
                    //cout<<"nstripsX="<<nStripInClusterX<<" | posX="<<posX<<" | X_MAXE="<<X_MAX_E_X<<" | minX="<<minPosX<<" | maxX="<<maxPosX<<" | mints="<<tsminx<<" | maxts="<<tsmaxx<<" | tswidth="<<(long long)tsmaxx-(long long)tsminx<<" | SumE = "<<E_X<<endl;
                    //! Y cluster here
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;

                    unsigned long long tsminy=0;
                    unsigned long long tsmaxy=0;
                    short minPosY=0;
                    short maxPosY=0;
                    double MAX_E_Y = -9999999;
                    short Y_MAX_E_Y = -1;
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        unsigned long long yts = yhitmaps_it->second->GetTimestamp();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y_prev;
                                //cout<<"nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        //!still in a same cluster
                        //if ( yen>0){
                            if (nStripInClusterY==0) {
                                minPosY=y;
                                tsminy=yts;
                            }
                            if (yen>MAX_E_Y){
                                MAX_E_Y=yen;
                                Y_MAX_E_Y = y;
                            }
                            if (yts<tsminy) tsminy=yts;
                            if (yts>tsmaxy) tsmaxy=yts;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                        //}
                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y;
                                //cout<<"last cluster: nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));


                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        y_prev = y;
                    }//y


                    nClusterX[z]++;
                }
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
                tsminx=0;
                tsmaxx=0;
                MAX_E_X = -9999999;
                X_MAX_E_X = -1;
            }
            //!still in a same cluster
            //if ( xen>0){
                if (nStripInClusterX==0) {
                    minPosX=x;
                    tsminx=xts;
                }
                if (xen>MAX_E_X){
                    MAX_E_X=xen;
                    X_MAX_E_X = x;
                }
                if (xts<tsminx) tsminx=xts;
                if (xts>tsmaxx) tsmaxx=xts;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
            //}
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    maxPosX = x;
                    //cout<<"last cluster: nstripsX="<<nStripInClusterX<<" | posX="<<posX<<" | X_MAXE="<<X_MAX_E_X<<" | minX="<<minPosX<<" | maxX="<<maxPosX<<" | mints="<<tsminx<<" | maxts="<<tsmaxx<<" | tswidth="<<(long long)tsmaxx-(long long)tsminx<<" | SumE = "<<E_X<<endl;
                    //! Y cluster here
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;

                    unsigned long long tsminy=0;
                    unsigned long long tsmaxy=0;
                    short minPosY=0;
                    short maxPosY=0;
                    double MAX_E_Y = -9999999;
                    short Y_MAX_E_Y = -1;
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        unsigned long long yts = yhitmaps_it->second->GetTimestamp();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y_prev;
                                //cout<<"nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */
                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        //!still in a same cluster
                        //if ( yen>0){
                            if (nStripInClusterY==0) {
                                minPosY=y;
                                tsminy=yts;
                            }
                            if (yen>MAX_E_Y){
                                MAX_E_Y=yen;
                                Y_MAX_E_Y = y;
                            }
                            if (yts<tsminy) tsminy=yts;
                            if (yts>tsmaxy) tsmaxy=yts;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                        //}
                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y;
                                //cout<<"last cluster: nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */
                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
                tsminx=0;
                tsmaxx=0;
                MAX_E_X = -9999999;
                X_MAX_E_X = -1;
            }
            x_prev = x;
        }//x

        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;
    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if (cluster->GetXTimestamp()>cluster->GetYTimestamp()){
                    cluster->SetTimestamp(cluster->GetYTimestamp());
                }else{
                    cluster->SetTimestamp(cluster->GetXTimestamp());
                }
                cluster->CalculateTimeDifference();

                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        ConstructEdiffRankingFlags();
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}

bool AIDA::IonGetPosAllNew2(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn

    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 10000000;
    //Int_t maxNposBeta = 10000000;
    Int_t maxNCluster = 10000000;

    Int_t nClusterX[NumDSSD];
    Int_t nClusterY[NumDSSD];
    memset(nClusterX,0,sizeof(nClusterX));
    memset(nClusterY,0,sizeof(nClusterY));

    //! Prepare maps of spatial distribution of hits
    //! NOTICE: remember to dallocate AIDAHit* !
    //! NOTICE: we use map: just select the first hit in each channel
    vector< map<short, AIDAHit* > > xhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator xhitmaps_it;
    vector< map<short, AIDAHit* > > yhitmaps(NumDSSD);
    map<short, AIDAHit* > ::iterator yhitmaps_it;

    for (size_t i=0;i<fhits.size();i++){
        AIDAHit* hit = fhits.at(i);
        int z = hit->GetZ();
        short xy = hit->GetXY();
        if (xy<128){
            xhitmaps[z].insert(std::make_pair(xy, hit));
        }else{
            yhitmaps[z].insert(std::make_pair(xy-128, hit));
        }
    }

    //cluster map sort by EX EY
    vector< multimap<Double_t,AIDACluster*> > cluster_map(NumDSSD);
    //multimap<E-corr, pair<Xhit,Yhit> >
    multimap<Double_t,AIDACluster*>::iterator cluster_map_it;

    Double_t posX=-11.;
    Double_t posY=-11.;
    //if (thereis) cout<<"----"<<endl;

    //! make them clusters!
    for (int z = 0;z < NumDSSD;z++){

        //! start with X clusters
        short x_prev = -11;
        short nStripInClusterX = 0;
        double E_X = 0;
        double E_X_ch = 0;
        unsigned long long tsminx=0;
        unsigned long long tsmaxx=0;
        short minPosX=0;
        short maxPosX=0;
        double MAX_E_X = -9999999;
        short X_MAX_E_X = -1;

        //! loop through the x axis
        for(xhitmaps_it = xhitmaps[z].begin(); xhitmaps_it != xhitmaps[z].end(); xhitmaps_it++){
            short x = xhitmaps_it->first;
            double xen = xhitmaps_it->second->GetEnergy();
            unsigned long long xts = xhitmaps_it->second->GetTimestamp();
            //!out of cluster
            if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    maxPosX = x_prev;
                    //cout<<"nstripsX="<<nStripInClusterX<<" | posX="<<posX<<" | X_MAXE="<<X_MAX_E_X<<" | minX="<<minPosX<<" | maxX="<<maxPosX<<" | mints="<<tsminx<<" | maxts="<<tsmaxx<<" | tswidth="<<(long long)tsmaxx-(long long)tsminx<<" | SumE = "<<E_X<<endl;
                    //! Y cluster here
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;

                    unsigned long long tsminy=0;
                    unsigned long long tsmaxy=0;
                    short minPosY=0;
                    short maxPosY=0;
                    double MAX_E_Y = -9999999;
                    short Y_MAX_E_Y = -1;
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        unsigned long long yts = yhitmaps_it->second->GetTimestamp();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y_prev;
                                //cout<<"nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        //!still in a same cluster
                        //if ( yen>0){
                            if (nStripInClusterY==0) {
                                minPosY=y;
                                tsminy=yts;
                            }
                            if (yen>MAX_E_Y){
                                MAX_E_Y=yen;
                                Y_MAX_E_Y = y;
                            }
                            if (yts<tsminy) tsminy=yts;
                            if (yts>tsmaxy) tsmaxy=yts;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                        //}
                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y;
                                //cout<<"last cluster: nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */

                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));


                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        y_prev = y;
                    }//y


                    nClusterX[z]++;
                }
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
                tsminx=0;
                tsmaxx=0;
                MAX_E_X = -9999999;
                X_MAX_E_X = -1;
            }
            //!still in a same cluster
            //if ( xen>0){
                if (nStripInClusterX==0) {
                    minPosX=x;
                    tsminx=xts;
                }
                if (xen>MAX_E_X){
                    MAX_E_X=xen;
                    X_MAX_E_X = x;
                }
                if (xts<tsminx) tsminx=xts;
                if (xts>tsmaxx) tsmaxx=xts;
                E_X += xen;
                E_X_ch += (double) x * xen;
                nStripInClusterX++;
            //}
            //! handle last cluster!
            if (xhitmaps_it == --xhitmaps[z].end() && nStripInClusterX > 0){
                if (nStripInClusterX<=maxStrInCluster){//copy to the last cluster handling in X
                    posX = E_X_ch/E_X;//center of gravity
                    maxPosX = x;
                    //cout<<"last cluster: nstripsX="<<nStripInClusterX<<" | posX="<<posX<<" | X_MAXE="<<X_MAX_E_X<<" | minX="<<minPosX<<" | maxX="<<maxPosX<<" | mints="<<tsminx<<" | maxts="<<tsmaxx<<" | tswidth="<<(long long)tsmaxx-(long long)tsminx<<" | SumE = "<<E_X<<endl;
                    //! Y cluster here
                    short y_prev = -11;
                    short nStripInClusterY = 0;
                    double E_Y = 0;
                    double E_Y_ch = 0;

                    unsigned long long tsminy=0;
                    unsigned long long tsmaxy=0;
                    short minPosY=0;
                    short maxPosY=0;
                    double MAX_E_Y = -9999999;
                    short Y_MAX_E_Y = -1;
                    //! loop through the y axis
                    for(yhitmaps_it = yhitmaps[z].begin(); yhitmaps_it != yhitmaps[z].end(); yhitmaps_it++){
                        short y = yhitmaps_it->first;
                        double yen = yhitmaps_it->second->GetEnergy();
                        unsigned long long yts = yhitmaps_it->second->GetTimestamp();
                        //!out of cluster
                        if (y - y_prev > 1 && y_prev != -11 && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y_prev;
                                //cout<<"nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */
                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        //!still in a same cluster
                        //if ( yen>0){
                            if (nStripInClusterY==0) {
                                minPosY=y;
                                tsminy=yts;
                            }
                            if (yen>MAX_E_Y){
                                MAX_E_Y=yen;
                                Y_MAX_E_Y = y;
                            }
                            if (yts<tsminy) tsminy=yts;
                            if (yts>tsmaxy) tsmaxy=yts;
                            E_Y += yen;
                            E_Y_ch += (double) y * yen;
                            nStripInClusterY++;
                        //}
                        //! handle last cluster!
                        if (yhitmaps_it == --yhitmaps[z].end() && nStripInClusterY > 0){
                            if (nStripInClusterY<=maxStrInCluster){//copy to the last cluster handling in X
                                posY = E_Y_ch/E_Y;//center of gravity
                                maxPosY = y;
                                //cout<<"last cluster: nstripsY="<<nStripInClusterY<<" | posY="<<posY<<" | Y_MAXE="<<Y_MAX_E_Y<<" | minY="<<minPosY<<" | maxY="<<maxPosY<<" | mints="<<tsminy<<" | maxts="<<tsmaxy<<" | tswidth="<<(long long)tsmaxy-(long long)tsminy<<" | SumE = "<<E_Y<<endl;
                                //! fill in data here
                                AIDACluster* cluster = new AIDACluster;
                                cluster->SetHitPosition(posX,posY,z);
                                cluster->SetXEnergy(E_X);
                                cluster->SetYEnergy(E_Y);
                                cluster->SetXMult(nStripInClusterX);
                                cluster->SetYMult(nStripInClusterY);
                                /*
                                //! Timming
                                short posXint = (short) round(posX);
                                short posYint = (short) round(posY);
                                AIDAHit* hitAtX = xhitmaps[z].find(posXint)->second;
                                AIDAHit* hitAtY = yhitmaps[z].find(posYint)->second;

                                //! Take the ealiest time stamp
                                if (hitAtX->GetFastTimestamp() > hitAtY->GetFastTimestamp()){
                                    cluster->SetTimestamp(hitAtY->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtY->GetFastTimestamp());
                                }else{
                                    cluster->SetTimestamp(hitAtX->GetTimestamp());
                                    cluster->SetFastTimestamp(hitAtX->GetFastTimestamp());
                                }
                                */
                                cluster->SetHitPositionMaxE(X_MAX_E_X,Y_MAX_E_Y,z);
                                cluster->SetXMinPos(minPosX);
                                cluster->SetXMaxPos(maxPosX);
                                cluster->SetYMinPos(minPosY);
                                cluster->SetYMaxPos(maxPosY);
                                cluster->SetXTimestamp(tsminx);
                                cluster->SetYTimestamp(tsminy);
                                cluster->SetXMaxTimestamp(tsmaxx);
                                cluster->SetYMaxTimestamp(tsmaxy);

                                //! finally insert map
                                if (E_X>E_Y) cluster_map[z].insert(make_pair(1-E_Y/E_X,cluster));
                                else cluster_map[z].insert(make_pair(1-E_X/E_Y,cluster));

                                nClusterY[z]++;
                            }
                            //! reset hits and energy in Y
                            E_Y = 0;
                            E_Y_ch = 0;
                            nStripInClusterY = 0;
                            tsminy=0;
                            tsmaxy=0;
                            MAX_E_Y = -9999999;
                            Y_MAX_E_Y = -1;
                        }
                        y_prev = y;
                    }//y

                    nClusterX[z]++;
                }
                //! reset hits and energy in X
                E_X = 0;
                E_X_ch = 0;
                nStripInClusterX = 0;
                tsminx=0;
                tsmaxx=0;
                MAX_E_X = -9999999;
                X_MAX_E_X = -1;
            }
            x_prev = x;
        }//x

        fnxclustersz[z]=nClusterX[z];
        if (nClusterX[z]>0) fnyclustersz[z]=nClusterY[z]/nClusterX[z];//nY cluster for each X clusters
    }//z

    Int_t maxZ=-1;
    for (int z = 0;z < NumDSSD;z++){
        //Record number of cluster here
        vector <Double_t> xindex;
        vector <Double_t> yindex;
        if (nClusterX[z]>0&&nClusterX[z]<maxNCluster&&nClusterY[z]>0&&nClusterY[z]<maxNCluster){
            for(cluster_map_it = cluster_map[z].begin(); cluster_map_it != cluster_map[z].end(); cluster_map_it++){
                AIDACluster* cluster = cluster_map_it->second;

                //! to fix bug on .9999999 number
                double xx=cluster->GetHitPositionX();
                double yy=cluster->GetHitPositionY();
                if(cluster->GetXMultiplicity()==1) xx=round(xx);
                if(cluster->GetYMultiplicity()==1) yy=round(yy);
                cluster->SetHitPosition(xx,yy,cluster->GetHitPositionZ());

                if (cluster->GetXTimestamp()>cluster->GetYTimestamp()){
                    cluster->SetTimestamp(cluster->GetYTimestamp());
                }else{
                    cluster->SetTimestamp(cluster->GetXTimestamp());
                }
                cluster->CalculateTimeDifference();

                if ((find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())){
                    cluster->SetRankingFlag(1);
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                }else{
                    if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut))) && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        delete cluster;
                    }
                }
                /*
                // old def.
                if (((corr_cut<=0)||((corr_cut>0)&&(cluster_map_it->first<corr_cut)))&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
                */
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        ConstructEdiffRankingFlags();
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}
