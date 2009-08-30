/*
gdc.h is library for transforming galaxy star catalog
to and from galaxy density data.

Copyright (C) 2009 Kresimir Cosic

gdc.h is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

You can contact the author by email  <cosic.kresimir@gmail.com>
*/



namespace gdc
{

typedef unsigned short uint16;
typedef unsigned char  uint8;

struct EncDecLogData
{
	uint16 maxDenKpc3Log;
	uint16 stepPercDiv;	

	double maxDenKpc3;
	double stepPerc;

	double densityFactor;

	int xBinN;
	int yBinN;
	int zBinN;
	float kpcBin;

	bool isOutsideRangeUp;
	bool isOutsideRangeDown;
	bool isOutsideRangeUpWarningDisplayed;
	bool isOutsideRangeDownWarningDisplayed;

	void ResetErrors(){
		isOutsideRangeUp=false;
		isOutsideRangeDown=false;
		isOutsideRangeUpWarningDisplayed=false;
		isOutsideRangeDownWarningDisplayed=false;
		}

	void SetSize(int xBinN,int yBinN,int zBinN, float kpcBin){
		this->xBinN=xBinN;
		this->yBinN=yBinN;
		this->zBinN=zBinN;
		this->kpcBin=kpcBin;
		}
	
	void SetParams(double totalDensityCalc, double starCount, double maxDenRelKpc3Requested,
			 double maxErrorRequested){

		ResetErrors();

		densityFactor=starCount/totalDensityCalc;

		double volume=double(xBinN)*yBinN*zBinN*kpcBin*kpcBin*kpcBin;
		double maxDenKpc3Requested=maxDenRelKpc3Requested*densityFactor;
		maxDenKpc3Log=uint16(10000+ ceil(log(maxDenKpc3Requested)/log(2.0)) );
		maxDenKpc3=pow(2.0,double(maxDenKpc3Log)-10000);

		#if VERBOSE_OUTPUT
			printf("maxDenRequested: %.0f  maxDen: %.0f\n",maxDenKpc3Requested,maxDenKpc3);
		#endif
		
		maxErrorRequested*=2;
		assert(maxErrorRequested>1.0/60000 && maxErrorRequested<=1.0);
		stepPercDiv=uint16( ceil(1/ maxErrorRequested ) );
		stepPerc=1.0/stepPercDiv;
		
		#if VERBOSE_OUTPUT
			printf("maxErrorRequested: %.5f  stepPerc: %.5f\n",maxErrorRequested/2,stepPerc);
		#endif				
		}	

	uint16 Encode(double densityRel) {
		double den=densityRel*densityFactor;
		double denBase=1+stepPerc;
		double halfDenBase=1+stepPerc/2;
		double val =32767 + log(den/maxDenKpc3) / log(denBase);
		
		if (val>32767*halfDenBase) {
			isOutsideRangeUp=true;
			if (!isOutsideRangeUpWarningDisplayed){
				isOutsideRangeUpWarningDisplayed=true;
				printf("Warning!!! : some data values outside up range, truncated to maxDensity\n");
				}
			}
		if (val<1/halfDenBase) {
			isOutsideRangeDown=true;
			if (!isOutsideRangeDownWarningDisplayed){
				isOutsideRangeDownWarningDisplayed=true;
				printf("Warning!!! : some data values outside down range, truncated to 0\n");
				}			
			}
		
		if (val<1/stepPerc) return 0;
		return uint16( floor(val+0.5) );		
		}
	double Decode(uint16 val) {		
		double denBase=1+stepPerc;						
		return maxDenKpc3*pow(denBase,int(val)-32767);		
		}
};

struct GdcFileWrite
{
	FILE* out;
	EncDecLogData par;
	
	bool isFirstVal;	

	bool shortDiffOn;
	std::vector<int> diffBuffer;	
	std::vector<int> range;	

	GdcFileWrite() {isFirstVal=true;} 	

	void WriteUint8(uint8 data)
	{				
		fwrite(&data,1,1,out);
		compressedSize+=1;
	}
	void WriteUint16(uint16 data)
	{
		WriteUint8(uint8(data>>8));
		WriteUint8(uint8(data&255));		
	}

	uint16 bitBuffer;
	int bitBufferSize;

	void WriteBits(int bits,int n){				
		uint16 mask=(1<<n)-1;
		uint16 val=uint16(bits)&mask;
		bitBuffer+=val<<bitBufferSize;
		bitBufferSize+=n;
		if (bitBufferSize>=8) {
			WriteUint8(bitBuffer&255);
			bitBuffer>>=8;
			bitBufferSize-=8;
			}	
		}
	void WriteFlushRemainingBits(){
		if (bitBufferSize>0) 
			WriteBits(0,8-bitBufferSize);
		}

	double totalSamples;
	double compressedSize;
	double statDiffBits[8];	
	int minVal,maxVal;

	void ResetStatistics() {
		totalSamples=0;
		compressedSize=0;
		for (int i=0;i<8;i++)
			statDiffBits[i]=0;
		minVal=0;maxVal=-10000000;
		}
	void PrintStatistics() {
		printf("  ");
		printf("TotSamples: %.0f  ",totalSamples);
		printf("OutBytes: %.0f  ",compressedSize);
		printf("maxVal:%6d  minVal:%6d\n",maxVal,minVal);
		#ifdef K_STAT
			printf("DiffBits:",compressedSize);
			for (int i=0;i<8;i++)
				printf(" %.0f",statDiffBits[i]);
			printf("\n");
		#endif
		}	
	
	int RoundBitsUp(int val,int bits){
		int mask=(1<<bits)-1;
		int valup=val+mask;
		return valup-(valup&mask);
		}

	int BitsPerSignedValue(int val)
	{
		if (val==0) return 0;
		for (int i=0;i<sizeof(int)*8-1;i++)
			if (val>=-(1<<i) && val<(1<<i) )
				return i+1;
					
		assert(false);
		return -1000;
	}

	void FlushDiffBuffer(){		
		while (diffBuffer.size()>0) {
			ProcessDiffBuffer();			
			}		
		}

	struct DynData
	{
		int bits;
		std::vector<int> len;	
		bool possible;
	};	
	std::vector<DynData> dynData;

	void DiffBufferBitDynamic()
	{
		dynData.resize(8);
		for (int i=0;i<8;i++) {
			dynData[i].bits=0;
			dynData[i].len.resize(1);
			dynData[i].len[0]=0;
			dynData[i].possible=true;
			}
		
		int n=int(range.size());		
		for (int i=0;i<n;i++) {
			int lev=range[i];
		
			int minBits=10000;
			int minJ=lev;

			for (int j=0;j<8;j++) {				
				if (!dynData[j].possible) continue;
								
				int bits=dynData[j].bits;
				if (bits<minBits) {
					minBits=bits;
					minJ=j;
					}					
				}
			std::vector<int>& minLen=dynData[minJ].len;
			
			for (int j=0;j<8;j++) {				
				if (j<lev) {
					dynData[j].possible=false;
					continue;
					}
				int bits1=dynData[j].bits+j;
				int bits2=minBits+9+j;

				if (bits1<bits2 && dynData[j].possible){
					dynData[j].bits=bits1;
					dynData[j].len.back()++;
					}
				else{
					dynData[j].bits=bits2;
					dynData[j].len=minLen;
					dynData[j].len.push_back(1);
					dynData[j].possible=true;
					}										
				}
			
			minJ++;
			}					
	}
	
	void WriteDiffBlock(int beg,int end){
		int len=end-beg;
		int bitRange=0;
		for (int j=beg;j<end;j++)
			bitRange=std::max(bitRange,range[j]);				
		
		assert(bitRange>=0 && bitRange<8);
		assert(len>=1 && len<=32);
		WriteBits(1,1);
		WriteBits(((len-1)<<3)+bitRange,8);
		
		for (int i=beg;i<end;i++)
			WriteBits(diffBuffer[i],bitRange);								
		}

	void ProcessDiffBuffer()
	{
		DiffBufferBitDynamic();

		int minI=0;
		int minBits=32*1000;
		for (int i=0;i<8;i++)
			if (dynData[i].possible && dynData[i].bits<minBits) {
				minI=i;
				minBits=dynData[i].bits;
				}

		std::vector<int>& parts=dynData[minI].len;				
		int partN=int(parts.size());
		if (partN>3) partN-=2;
		else if (partN>1) partN=1;		

		int beg=0;
		for (int i=0;i<partN;i++){
			int end=beg+parts[i];			
			WriteDiffBlock(beg,end);
			beg=end;
			}

		int len=beg;
		int diffBufferN=int(diffBuffer.size());
		for (int i=0;i<diffBufferN-len;i++) {
			diffBuffer[i]=diffBuffer[i+len];
			range[i]     =range[i+len];
			}
		diffBuffer.resize(diffBufferN-len);
		range.resize(diffBufferN-len);
	}	

	void WriteGdcData(double denRel)
	{		
		static uint16 val=0;
		static int rleDiff=0;
		
		uint16 oldVal=val;
		val=par.Encode(denRel);
		totalSamples++;
		if (maxVal<val-32767) maxVal=val-32767;
		if (minVal>val-32767) minVal=val-32767;

		if (isFirstVal) {
			isFirstVal=false;		
			WriteUint16(val);			
			shortDiffOn=false;			
			rleDiff=0;
			bitBufferSize=0;
			bitBuffer=0;
			return;
			}
		
		int oldRleDiff=rleDiff;
		rleDiff=(signed short int)val-oldVal;

		#ifdef K_STAT
			int rleDiffError=rleDiff-oldRleDiff;
			int bitsRleDiffError=BitsPerSignedValue(rleDiffError);
			if (bitsRleDiffError>7) bitsRleDiffError=7;
			statDiffBits[bitsRleDiffError]++;
		#endif

		if (!shortDiffOn){			
			if (rleDiff>=-64 && rleDiff<64) {
				WriteUint8( rleDiff|128 );
				shortDiffOn=true;				
				diffBuffer.resize(0);
				}
			else
				WriteUint16(val);
						
			return;			
			}
						
		int rleD2=rleDiff-oldRleDiff;
		if ( !(rleD2>=-64 && rleD2<64) ) {
			FlushDiffBuffer();
			WriteBits(0,1);
			WriteFlushRemainingBits();
			WriteUint16(val);			
			shortDiffOn=false;			
			return;
			}
		
		diffBuffer.push_back(rleD2);
		range.push_back(BitsPerSignedValue(rleD2));
		if (diffBuffer.size()>=32)
			ProcessDiffBuffer();	
	}
	void WriteGdcDataFlush(){
		if (shortDiffOn) {
			FlushDiffBuffer();
			WriteFlushRemainingBits();
			}
		}
	void WriteGdcHeader()
	{
		WriteUint8('G'); //readable file marker
		WriteUint8('D');
		WriteUint8('C');
		WriteUint8(' ');
		WriteUint8('F');
		WriteUint8('I');
		WriteUint8('L');
		WriteUint8('E');	
		WriteUint8(13);
		WriteUint8(10);

		WriteUint8(1);   //major version
		WriteUint8(0);   //minor version
		WriteUint16(par.xBinN);
		WriteUint16(par.yBinN);
		WriteUint16(par.zBinN);
		
		uint16 xSizeKpc=uint16(par.xBinN*par.kpcBin+0.5);
		if (xSizeKpc!=par.xBinN*par.kpcBin) {
			printf("Warning: could not write exact size of galaxy!\n");
			printf("Requested size: %.4f   Approximated by: %d!\n",par.xBinN*par.kpcBin,xSizeKpc);
			}
		
		WriteUint16(xSizeKpc);
		WriteUint16(par.maxDenKpc3Log);
		WriteUint16(par.stepPercDiv);	
	}

}; //GdcFileWrite ends

struct GdcFileRead
{
	FILE* in;
	EncDecLogData par;
	
	bool isFirstVal;		

	GdcFileRead() {isFirstVal=true;} 	

	size_t bytesRead;
	double totalDensity;

	uint16 ReadUint16()
	{
		uint16 ret;
		ret=uint16(ReadUint8())<<8;
		ret+=ReadUint8();
		return ret;
	}
private:
	uint8 ReadUint8()
	{				
		uint8 ret;
		fread(&ret,1,1,in);
		if (feof(in)) {
			printf("\nUnexpected end of file!\n");
			abort();
			}
		bytesRead++;
		return ret;
	}
	

	uint16 bitBuffer;
	int bitBufferSize;
	
	int ReadBits(int n){				
		assert(n>=0 && n<=8);
		
		uint16 mask=(1<<n)-1;

		if (bitBufferSize<n) {						
			uint16 data=ReadUint8();
			data<<=bitBufferSize;
			bitBuffer+=data;
			bitBufferSize+=8;							
			}	
		bitBufferSize-=n;
		uint16 ret=bitBuffer&mask;
		bitBuffer>>=n;

		if (n>0) {
			bool isSigned=(ret&(1<<(n-1)))!=0;
			if (isSigned)
				return int(ret)-(1<<n);
			}

		return ret;
		}	
	
	struct GdcDataState
	{		
		int rleDiff;
		bool shortDiffOn;
		int bitRange;
		int bitRangeLen;
		uint16 val;

		void SetRleDiff(int p_rleDiff){
			rleDiff=p_rleDiff;
			int bigVal=val;
			bigVal+=rleDiff;
			if ((bigVal&32767) != bigVal){
				printf("\nData error: packed value outside range!\n");
				abort();
				}
			val=uint16(bigVal);			
			}
	};

	uint16 GdcReadPacked(GdcDataState& gds)
	{
		int packed=ReadBits(gds.bitRange);		
		gds.bitRangeLen--;

		gds.SetRleDiff(gds.rleDiff+packed);		
		return gds.val;
	}
	void GdcReadTag(GdcDataState& gds)
	{		
		uint8 tag=uint8(ReadBits(8));
		gds.bitRangeLen=(tag>>3)+1;		
		gds.bitRange=tag&7;	
	}

	uint16 ReadGdcDataVal()
	{		
		static GdcDataState gds;
		
		if (!isFirstVal && gds.shortDiffOn) {
			if (gds.bitRangeLen!=0) 
				return GdcReadPacked(gds);
			
			bool isShort=ReadBits(1)!=0;
			if (isShort) {
				GdcReadTag(gds);
				return GdcReadPacked(gds);
				}	
			else
				{
				gds.shortDiffOn=false;
				bitBufferSize=0;
				bitBuffer=0;
				}
			}

		uint8 byte1=ReadUint8();
		bool isRelative=(byte1&128)!=0;		

		if (!isFirstVal && !gds.shortDiffOn){			
			if (!isRelative) {
				uint8 byte2=ReadUint8();
				gds.val=(uint16(byte1)<<8)+byte2;
				return gds.val;
				}
			
			gds.shortDiffOn=true;						
			gds.SetRleDiff(((signed char)(byte1<<1))>>1);		
			return gds.val;			
			}
				
		uint8 byte2=ReadUint8();
		gds.val=(uint16(byte1)<<8)+byte2;		

		assert(isFirstVal);			
		if (isRelative){
			printf("\nData error: first value can not be relative!\n");
			abort();
			}
		isFirstVal=false;
		totalDensity=0;
		gds.shortDiffOn=false;
		gds.bitRangeLen=0;
		gds.bitRange=0;
		gds.rleDiff=0;
		bitBufferSize=0;
		bitBuffer=0;
		return gds.val;
	}


public:
	uint16 minVal,maxVal;		
	
	void ResetStatistics() {		
		minVal=-1;maxVal=0;
		}
	void PrintStatistics() {
		double maxDensity=par.Decode(maxVal);
		double minDensity=par.Decode(minVal);
		printf("maxVal:%6d  minVal:%6d  maxDen:%11.0f  minDen: %.10f\n",
					int(maxVal)-32767,int(minVal)-32767,maxDensity, minDensity);		
		}	

	#if TEST_CHECK_GDC_COMPRESSION
		uint16 prevVal;
	#endif
	double ReadGdcData()
	{
		uint16 val=ReadGdcDataVal();
		if (maxVal<val) maxVal=val;
		if (minVal>val) minVal=val;		

		#if TEST_CHECK_GDC_COMPRESSION
			prevVal=val;;
		#endif
		double den=par.Decode(val);
		totalDensity+=den;
		return den;
	}

	int ReadGdcHeader()
	{
		if (ReadUint8()!='G') return 1;
		if (ReadUint8()!='D') return 2;
		if (ReadUint8()!='C') return 3;
		if (ReadUint8()!=' ') return 4;
		if (ReadUint8()!='F') return 5;
		if (ReadUint8()!='I') return 6;
		if (ReadUint8()!='L') return 7;
		if (ReadUint8()!='E') return 8;
		if (ReadUint8()!=13) return 9;
		if (ReadUint8()!=10) return 10;						

		if (ReadUint8()!=1) return 11; //major version
		if (ReadUint8()!=0) return 12; //minor version

		par.xBinN=ReadUint16();
		par.yBinN=ReadUint16();
		par.zBinN=ReadUint16();				
		
		par.kpcBin=float(ReadUint16())/par.xBinN;
		par.maxDenKpc3Log=ReadUint16();
		par.stepPercDiv=ReadUint16();

		par.maxDenKpc3=pow(2.0,double(par.maxDenKpc3Log)-10000);		
		par.stepPerc=1.0/par.stepPercDiv;

		if (feof(in)) return 20;

		return 0;								
	}

}; //GdcFileRead ends

} //namespace gdc ends

