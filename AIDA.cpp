#include "AIDA.h"

bool AIDA::BetaGetPosNew(Double_t corr_cut,Double_t sumexcut[],Double_t sumeycut[])
{
    //!CAUTION: condition if (x - x_prev > 1 && x_prev != -11 && nStripInClusterX > 0){
    //! is apply so we should assume "good" calibrationn


    //! A  tiny event builder by dimension!
    //! So far the method can only handle the first hit in each strip!
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

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
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
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

                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;
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
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

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
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
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



                /*
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;

                }else{
                    //! Deallocating memory
                    delete cluster;
                }
                */
                if ((corr_cut<=0)||(corr_cut>0)&&(cluster_map_it->first<corr_cut)&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
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
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

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
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
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

                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;
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
    Int_t maxStrInCluster = 1000;
    Int_t maxNposBeta = 1000;
    Int_t maxNCluster = 1000;

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
    }//z

    Int_t maxZ=-1;

    for (int z = 0;z < NumDSSD;z++){
        Int_t mult_z = 0;
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



                /*
                if (((corr_cut<=0)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end()))||
                        ((corr_cut>0)&&(cluster_map_it->first<corr_cut)&&(find(xindex.begin(),xindex.end(),cluster->GetHitPositionX())==xindex.end())&&(find(yindex.begin(),yindex.end(),cluster->GetHitPositionY())==yindex.end())))
                {
                    //cout<<cluster_map_it->first<<"-"<<cluster->GetXEnergy()<<"-"<<cluster->GetYEnergy()<<endl;
                    if (mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                        this->AddCluster(cluster);
                        maxZ = z;
                    }else{
                        //! Deallocating memory
                        delete cluster;
                    }
                    xindex.push_back(cluster->GetHitPositionX());
                    yindex.push_back(cluster->GetHitPositionY());
                    mult_z++;

                }else{
                    //! Deallocating memory
                    delete cluster;
                }
                */
                if ((corr_cut<=0)||(corr_cut>0)&&(cluster_map_it->first<corr_cut)&&mult_z<maxNposBeta && cluster->GetXEnergy()>sumexcut[z] && cluster->GetYEnergy()>sumeycut[z]){
                    this->AddCluster(cluster);
                    maxZ = z;
                }
                else delete cluster;
            }//loop all cluster map
        }//if condition
    }//all z


    if(this->GetNClusters()>0) {
        this->SetMaxZ(maxZ);
        return true;
    }
    return false;
}
