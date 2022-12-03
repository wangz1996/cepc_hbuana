#include "DatManager.h"
#include "Global.h"
using namespace std;
extern char char_tmp[200];
extern int int_tmp=0;
int DatManager::CatchSPIROCBag(ifstream &f_in, vector<int> &buffer_v, int &layer_id,int &cycleID,int &triggerID){
    //cout<<"catch a bag"<<endl;
    bool b_begin=0;
    bool b_end=0;
    int buffer=0;
    while(!b_begin && f_in.read((char*)(&buffer),1) ){
        //cout<<hex<<buffer<<" ";
        buffer_v.push_back(buffer);
        if(buffer_v.size()>4) buffer_v.erase(buffer_v.begin(),buffer_v.begin()+buffer_v.size()-4);
        if(buffer_v[0]==0xfa && buffer_v[1]==0x5a && buffer_v[2]==0xfa && buffer_v[3]==0x5a && buffer_v.size()==4) b_begin=1;
    }
    while(!b_end && f_in.read((char*)(&buffer),1)){
        buffer_v.push_back(buffer);
        int_tmp = buffer_v.size();
        //if(int_tmp>4 && buffer_v[int_tmp-2] == 0xfe && buffer_v[int_tmp-1] == 0xee && buffer_v[int_tmp-4] == 0xfe && buffer_v[int_tmp-3] == 0xee) b_end=1;
        if(int_tmp>=4 && buffer_v[int_tmp-2] == 0xfe && buffer_v[int_tmp-1] == 0xee && buffer_v[int_tmp-4] == 0xfe && buffer_v[int_tmp-3] == 0xee) b_end=1;
    }
    f_in.read((char*)(&buffer),1); 
    if(buffer!=0xff){
        cout<<" abnormal layer ff "<<hex<<buffer<<endl;
        buffer_v.clear();
        return 0;
    }
    f_in.read((char*)(&buffer),1); 
    if(buffer<0 || buffer>39){
        cout<<" abnormal layer "<<hex<<buffer<<endl;
        buffer_v.clear();
        return 0;
    }
    layer_id=buffer;        
    //cout<<"cycleID "<<hex<<cycleID<<endl;
    if( (buffer_v.size())%2 ){
        cout<<"wrong bag size "<<dec<<buffer_v.size()<<endl;
        buffer_v.clear();
        return 0;//b_readover
    }
    for (int i = 0; i < buffer_v.size()/2; ++i){
        buffer_v[i] = buffer_v[2*i]*0x100+buffer_v[2*i+1];
    }
    buffer_v.erase(buffer_v.begin()+buffer_v.size()/2,buffer_v.end());
    cycleID = buffer_v[2]*0x10000+buffer_v[3];
    triggerID = buffer_v[4];
    buffer_v.erase(buffer_v.begin()+2,buffer_v.begin()+5);
    if(f_in.eof())return 0;
    else return 1;
}
int DatManager::DecodeAEvent(vector<int> &chip_v,int layer_id,int Memo_ID){
    int size=chip_v.size();
    if((size%73)!=1){
        cout<<"wrong chip size "<<chip_v.size()<<endl;
        return 0;
    }
    int Memo_No=size/channel_FEE;
    int offset=72*(Memo_No-1);    
    for(int i_ch=0; i_ch<channel_No; ++i_ch){//72+BCID+Chip
        int chanID = channel_No-1-i_ch;
        int gain = chip_v.at(i_ch+offset)&0x2000;
        int hit  = chip_v.at(i_ch+offset)&0x1000;
        int tdc  = chip_v.at(i_ch+offset)&0x0fff;
        int adc  = chip_v.at(i_ch+offset+channel_No)&0x0fff;
        int chip = chip_v[size-1]-1;
        int BCID = chip_v[size-2];           
        //if(hit<=1)continue; 
        _cellID.push_back(layer_id*1E5+chip*pow(10,4)+Memo_ID*pow(10,2)+chanID);
        _bcid.push_back(BCID);
        if(hit>1)_hitTag.push_back(1);
        else _hitTag.push_back(0);

        if(gain>1)_gainTag.push_back(1);
        else _gainTag.push_back(0);
        _charge.push_back(adc);
        _time.push_back(tdc);
    }
    /*for (int i = 0; i < chip_v.size(); ++i){
        cout<<hex<<chip_v[i]<<" ";
    }
    cout<<endl;
    */
    chip_v.erase(chip_v.end()-2,chip_v.end()-1);
    chip_v.erase(chip_v.begin()+offset,chip_v.begin()+offset+72);
    if(chip_v.size()==1)chip_v.clear();
}

int DatManager::FillChipBuffer(vector<int> &buffer_v,int cycleID,int triggerID,int layer_id){
    int size = buffer_v.size();
    if(size<4){
        //cout<<"FillChipBuffer:wrong bag size "<<size<<endl;
        buffer_v.clear();
        return 0;
    }
    if( buffer_v[0]!=0xfa5a || buffer_v[1]!=0xfa5a || buffer_v[size-2]!=0xfeee || buffer_v[size-1]!=0xfeee){
        cout<<"FillChipBuffer:wrong bag package "<<hex<<buffer_v[0]<<" "<<buffer_v[1]<<" "<<buffer_v[size-2]<<" "<<buffer_v[size-1]<<endl;
        buffer_v.clear();
        return 0;
    }
    else{
        buffer_v.erase(buffer_v.begin(),buffer_v.begin()+2);
        buffer_v.erase(buffer_v.end()-2,buffer_v.end());
        size = buffer_v.size();
    }
    for (int i = 0; i < buffer_v.size(); ++i){
        //cout<<hex<<buffer_v[i]<<" ";
    }
    for (int i = channel_FEE; i<buffer_v.size(); i=i+channel_FEE){
        //cout<<dec<<i<<" "<<buffer_v.size()<<" "<<hex<<buffer_v[i]<<endl;
        if (buffer_v[i]<1 || buffer_v[i]>9) continue;
        int chip=buffer_v[i]-1;
        _chip_v[layer_id][chip].assign(buffer_v.begin(),buffer_v.begin()+i+1);
        buffer_v.erase(buffer_v.begin(),buffer_v.begin()+i+1);
        //cout<<endl<<dec<<_chip_v[layer_id][chip].back()<<" FillChipBuffer "<<" "<<buffer_v.size()<<endl;
        i=0;
    }
    if(buffer_v.size()){
        count_chipbuffer++;
        cout<<hex<<cycleID<<" "<<buffer_v.back()<<" FillChipBuffer:abnormal chip buffer "<<dec<<" "<<layer_id<<" "<<buffer_v.size()<<" "<<count_chipbuffer<<endl;
        for (int i_chip = 0; i_chip < chip_No; ++i_chip){
            //_chip_v[layer_id][i_chip].clear();
        }
        buffer_v.clear();
        return 0;
    }
    return 1;
}

int DatManager::Decode(const string& input_file,const string& output_file)
{
    ifstream f_in;
    int layer_id;
    int cycleID;
    int triggerID;
    int BCID[Layer_No][chip_No];
    int Memo_ID[Layer_No][chip_No];
    f_in.open(input_file,ios::in);
    if(!f_in){
        cout<<"cant open "<<input_file<<endl;
        return 0;
    }
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer){
        _buffer_v.clear();
        for (int i_chip = 0; i_chip < chip_No; ++i_chip){
            BCID[i_layer][i_chip]=-1;
            Memo_ID[i_layer][i_chip]=0;
            _chip_v[i_layer][i_chip].clear();
        }
    }
    //string str_out=outputDir+"/"+"cosmic.root";
    string tmp_outname=input_file;
	tmp_outname=tmp_outname.substr(tmp_outname.find_last_of('/')+1);
	tmp_outname=tmp_outname.substr(0,tmp_outname.find_last_of('.'));
    string str_out=output_file+"/"+tmp_outname+".root";
	std::cout<<str_out.c_str()<<std::endl;
    TFile *fout;
    fout =  TFile::Open(str_out.c_str(),"RECREATE");
    if(!fout){
        cout<<"cant create "<<str_out<<endl;
        return 0;
    }
    TTree *tree = new TTree("Cosmic_Event","data from binary file");
    SetTreeBranch(tree);
    int Bag_No=0;
    int Event_No=0;
    int Abnormal_Event_No=0;
    int coincidence_No=0;
    long pre_trigID=0;
    long pre_cycleID=0;
    long tmp_ID=0;
    bool b_ReadOver=1;
    bool b_chipbuffer=0;
    while(!(f_in.eof()) || b_chipbuffer){
    //while((!(f_in.eof()) || b_chipbuffer) && Event_No<=1E2){
        _buffer_v.clear();
        CatchSPIROCBag(f_in,_buffer_v,layer_id,cycleID,triggerID);
        b_chipbuffer=Chipbuffer_empty();
        if(_buffer_v.size()<74){
            if(_buffer_v.size()!=4)cout<<"abnormal SPIROC bag size "<<_buffer_v.size()<<endl;
            _buffer_v.clear();
            if(!b_chipbuffer) continue;
        }
        //tmp_ID=long(cycleID);
        tmp_ID=long(triggerID);               
        if(!b_chipbuffer) {
            pre_trigID=triggerID;
            pre_cycleID=cycleID;
        }
        if(triggerID<pre_trigID  && fabs( triggerID - pre_trigID ) <10000){
            cout<<pre_cycleID<<" "<<pre_trigID<<" abnormal ID "<<cycleID<<" "<<triggerID<<endl;
            _buffer_v.clear();
        }
        //cout<<hex<<triggerID<<" tmp_ID "<<tmp_ID<<endl;
        Bag_No++;
        if((Bag_No%10000)==0)cout<<" Bag_No "<<Bag_No<<endl;
        if(triggerID==pre_trigID){
            FillChipBuffer(_buffer_v,cycleID,triggerID,layer_id);
        }    
        if (triggerID>pre_trigID || f_in.eof() || fabs( triggerID - pre_trigID ) >40000){
            if( fabs( triggerID - pre_trigID ) >40000)cout<<hex<<triggerID<<" Debug "<<pre_trigID<<endl;
            if( ( triggerID - pre_trigID ) !=1 ) Abnormal_Event_No++;
            while(b_chipbuffer){
                long EventID=0x1000000000000;
                int b_layer=0;
                for (int i_layer = 0; i_layer < Layer_No; ++i_layer)
                {
                    for (int i_chip = 0; i_chip < chip_No; ++i_chip){
                        int size=_chip_v[i_layer][i_chip].size();
                        if(size==0)continue;
                        //cout<<hex<<size<<" "<<" "<<i_chip<<" "<<i_layer<<endl;
                        Memo_ID[i_layer][i_chip]=size/73;
                        BCID[i_layer][i_chip]=_chip_v[i_layer][i_chip].at(size-2);
                        tmp_ID=long(BCID[i_layer][i_chip]);
                        if(EventID>tmp_ID){
                            EventID=tmp_ID;
                        }
                    }
                }
                int b_coincidence=0;
                BranchClear();
                for (int i_layer = 0; i_layer < Layer_No; ++i_layer){
                    for (int i_chip = 0; i_chip < chip_No; ++i_chip){
                        int size=_chip_v[i_layer][i_chip].size();
                        long ID_tmp=long(BCID[i_layer][i_chip]);
                        //cout<<hex<<size<<" "<<ID_tmp<<" "<<i_chip<<" "<<i_layer<<endl;
                        if(size==0)continue;
                        /*if(EventID!=ID_tmp){
                            cout<<ID_tmp<<" different ID "<<EventID<<" "<<i_layer<<" "<<i_chip<<endl;
                        }
                        else cout<<ID_tmp<<" same ID "<<EventID<<" "<<i_layer<<" "<<i_chip<<endl;*/
                        if(1){    
                            //cout<<hex<<ID_tmp<<" coincidence "<<EventID<<dec<<" "<<i_chip<<" "<<i_layer<<endl;
                            b_coincidence++;
                            _cycleID=pre_cycleID;
                            _triggerID=pre_trigID;
                            DecodeAEvent(_chip_v[i_layer][i_chip],i_layer,Memo_ID[i_layer][i_chip]-1);
                            Memo_ID[i_layer][i_chip]--;
                            _chip_v[i_layer][i_chip].clear();
                        }
                    }
                }
                Event_No++;
                //if(b_coincidence==2){coincidence_No++;cout<<"coincidence No "<<coincidence_No<<endl;}
                //if((_cycleID-pre_trigID)!=1) {coincidence_No++;cout<<hex<<_cycleID<<" cycle error "<<pre_trigID<<endl;}
                tree->Fill();
                BranchClear();
                b_chipbuffer=Chipbuffer_empty();
                if(b_chipbuffer)cout<<b_chipbuffer<<" Chip buffer not clear "<<pre_trigID<<endl;
                //cout<<Event_No<<" Event "<<pre_trigID<<" "<<b_chipbuffer<<endl;
            }
            //cout<<b_chipbuffer<<" Chip buffer clear "<<pre_trigID<<" "<<triggerID<<endl;
            FillChipBuffer(_buffer_v,cycleID,triggerID,layer_id);
            pre_trigID=triggerID;
            pre_cycleID=cycleID;                
        }    
    }
    cout<<dec<<Abnormal_Event_No<<" Abnormal Event No in "<<Event_No<<endl;
    f_in.close();
    tree->Write();
    fout->Write();
    fout->Close();
    return 1;
}

void DatManager::SetTreeBranch(TTree *tree){
	tree ->Branch("cycleID",&_cycleID);
	tree ->Branch("triggerID",&_triggerID);
	tree ->Branch("cellIDs",&_cellID);
	tree ->Branch("BCIDs",&_bcid);
	tree ->Branch("hitTags",&_hitTag);
	tree ->Branch("gainTags",&_gainTag);
	tree ->Branch("charges",&_charge);
	tree ->Branch("times",&_time);
}

void DatManager::BranchClear() 
{
	_cellID.clear();
	_bcid.clear();
	_hitTag.clear();
	_gainTag.clear();
	_charge.clear();
	_time.clear();    
}
/*int raw2Root::RMFelixTag(string inputDir,string outputDir){
    ifstream f_datalist,f_in[Layer_No];
    ofstream f_out[Layer_No];
    string str_datalist=inputDir+"/datalist";
    f_datalist.open(str_datalist,ios::in);
    if (!f_datalist){
        cout<<"cant open "<<str_datalist<<endl;
        return 0;
    }
    for (int i = 0; f_datalist>>str_tmp; ++i){
        f_in[i].open(str_tmp,ios::in);
        if (!f_in[i]){
            cout<<"f_in cant open "<<str_tmp<<endl;
            continue;
        }
        else cout<<"Read "<<str_tmp<<endl;
        str_tmp=find_datname(str_tmp);
        str_tmp=outputDir+"/"+str_tmp+".dat";
        cout<<str_tmp<<endl;
        f_out[i].open(str_tmp,ios::out|ios::binary);
        if (!f_out[i]){
            cout<<"f_out cant open "<<str_tmp<<endl;
            continue;
        }   
        else cout<<"Write "<<str_tmp<<endl;
        bool b_felix=0;
        int buffer=0;
        int size=0;
        int count=0;
        std::vector<int> v;
        //ffff 0000 baba 5afff ******* ffa5 abab 0000 ffff
        while(f_in[i].read((char*)(&buffer),1)){
            count++;
            v.push_back(buffer);
            int size = v.size();
            if( count>=6 && count%1024 == 4 && v[size-1] == 0xab && v[size-2] == 0xcd){
                v.erase(v.end()-6,v.end());
                continue;
            }
            if( size%2==1 || size<8 )continue;
            if( v[size-1]==0xff && v[size-2]==0xff && v[size-3]==0x00 && v[size-4]==0x00 &&
                v[size-5]==0xab && v[size-6]==0xab && v[size-7]==0xa5 && v[size-8]==0xff){//ffa5 abab 0000 ffff
                v.clear();
                v.push_back(0x00);v.push_back(0x00);v.push_back(0xff);v.push_back(0xff);//0000 ffff
                continue;
            }
            if( v[size-1]==0xff && v[size-2]==0x5a && v[size-3]==0xba && v[size-4]==0xba &&
                v[size-5]==0x00 && v[size-6]==0x00 && v[size-7]==0xff && v[size-8]==0xff){//ffff 0000 baba 5aff
                if(v[0]==0x00 && v[1]==0x00 && v[2]==0xff && v[3]==0xff) v.erase(v.begin(),v.begin()+4);
                v.erase(v.end()-8,v.end());
                size=v.size();
                /*if(v[size-20]==0xfa && v[size-19]==0x23 && v[size-16]==0xcd && v[size-15]==0xab){
                    v.erase(v.end()-20,v.end()-14);
                    size=v.size();
                }
                else cout<<"no felix tag before"<<hex<<v[size-20]<<" "<<v[size-19]<<" "<<v[size-18]<<" "<<v[size-17]<<" "<<v[size-16]<<" "<<endl;
                for (int i_v = 0; i_v < size; ++i_v){
                    buffer=v[i_v];
                    f_out[i].write((char*)(&buffer),1);
                }
                v.clear();                
                v.push_back(0xff);v.push_back(0xff);v.push_back(0x00);v.push_back(0x00);//ffff 0000
                continue;
            }
        }
        if(v[0]==0xff && v[1]==0xff && v[2]==0x00 && v[3]==0x00) v.clear();
        else if(v[0]==0x00 && v[1]==0x00 && v[2]==0xff && v[3]==0xff){
            for (int i_v = 4; i_v < v.size(); ++i_v){
                buffer=v[i_v];
                f_out[i].write((char*)(&buffer),1);
            }            
        }
        else cout<<"Felix Tag error "<<hex<<v[0]<<v[1]<<v[2]<<v[3]<<endl;
        v.clear();
    } 
    return 1;
}*/
DatManager::~DatManager()
{}
