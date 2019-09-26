@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_WATER_1
//#define USE_CLOUDS_1
@import ./scene;

// textureMaps[0] = 'images/raymarch/organic2.jpg' (1024 x 1024)
// textureMaps[1] = 'images/raymarch/abstract1.jpg' (1024 x 1024)
// textureMaps[2] = 'images/raymarch/grayNoise256.png' (256 x 256)
// textureMaps[3] = 'images/raymarch/lichen.jpg' (1024 x 1024)
#define iChannel0   textureMaps[0]
#define iChannel1   textureMaps[1]
#define iChannel2   textureMaps[2]
#define iChannel3   textureMaps[3]

// Choose...
#define ENABLE_AUTO_VIEW
#define ENABLE_FISH_COLORING

const float FAR = 300.0;
const vec3 SUN_LIGHT = normalize(vec3(0.9,0.35,-0.2));
const vec3 SUN_COLOR = vec3(1.0, 0.9, 0.85);
const float SKY_HEIGHT = 100.0;
const vec3 WATER_COLOR = vec3(0.3, 0.4, 0.45);
const vec3 UNDER_WATER_COLOR = vec3(0.0, 0.15, 0.25);

float skyMap( in vec3 p )
{
    return SKY_HEIGHT - p.y;
}

float aquaTerrainH( vec3 p )
{
    float h = 1.0;
	vec3 q = p;
	float th = smoothstep( 0.1, 0.4, textureLod( iChannel0, 0.002*q.xz, 0.0 ).x );
    float rr = smoothstep( 0.2, 0.5, textureLod( iChannel1, 2.0*0.02*q.xz, 0.0 ).y );
	h = 0.9 + (1.0-0.6*rr)*(1.5-1.0*th) * 0.1*(1.0-textureLod( iChannel0, 0.1*q.xz, 0.0 ).x);
	h += th * 1.25;
    h -= 0.24 * rr;
	h *= 0.75;
    return -h;
}

// float aquaTerrainL( in sampler2D tex, in vec3 p )
// {
//     return aquaTerrainH( p );
// }

float aquaTerrainMap( vec3 p )
{
    return (p.y - aquaTerrainH(p))*1.0;//0.3;
}

vec3 aquaTerrainColor( in vec3 rd, in vec3 pos, in vec3 nor )
{
    vec4 mate = vec4(0.5,0.5,0.5,0.0);

    vec3 te = texture( iChannel0, 0.1*pos.xz ).xyz;
    te = 0.05 + te;
    mate.xyz = 0.6*te;          // material diffuse
    mate.w = 5.0*(0.5+0.5*te.x);// material shininess

    float th = smoothstep( 0.1, 0.4, texture( iChannel0, 0.002*pos.xz ).x );
    vec3 dcol = mix( vec3(0.1, 0.1, 0.0), 0.4*vec3(0.65, 0.4, 0.2), 0.2+0.8*th );
    mate.xyz = mix( mate.xyz*0.5, dcol, th*smoothstep( 0.0, 1.0, nor.y ) );

    float rr = smoothstep( 0.2, 0.4, texture( iChannel1, 2.0*0.02*pos.xz ).y );
    mate.xyz *= mix( vec3(1.0), vec3(0.2,0.2,0.2)*1.5, rr );
    mate.xyz *= 1.5;

    // lighting
    vec3 lig = SUN_LIGHT;
    float sky = saturate(nor.y);
    float bou = saturate(-nor.y);
    float dif = max(dot(nor,lig),0.0);
    float bac = max(0.3 + 0.7*dot(nor,-vec3(lig.x,0.0,lig.z)),0.0);
    float sha = 0.0; if( dif>0.001 ) sha = sceneShadow( pos+0.01*nor, lig, 0.0005, 32.0 );
    float fre = pow( saturate( 1.0 + dot(nor,rd) ), 5.0 );
    float spe = max( 0.0, pow( saturate( dot(lig, reflect(rd,nor)) ), mate.w ) ) * mate.w;
    float sss = pow( saturate( 1.0 + dot(nor,rd) ), 3.0 );
    
    // water caustics...
    float cc  = 0.55*texture( iChannel2, 1.8*0.02*pos.xz + 0.007*time*vec2( 1.0, 0.0) ).x;
          cc += 0.25*texture( iChannel2, 1.8*0.04*pos.xz + 0.011*time*vec2( 0.0, 1.0) ).x;
          cc += 0.10*texture( iChannel2, 1.8*0.08*pos.xz + 0.014*time*vec2(-1.0,-1.0) ).x;
    cc = 0.6*(1.0-smoothstep( 0.0, 0.025, abs(cc-0.4))) + 0.4*(1.0-smoothstep( 0.0, 0.150, abs(cc-0.4)));
    dif *= 1.0 + 2.0*cc;

    // lights
    vec3 lin = vec3(0.0);
    lin += 3.5*dif*vec3(1.00,1.00,1.00)*sha;
    lin += 3.0*sky*vec3(0.10,0.20,0.35);
    lin += 1.0*bou*vec3(0.20,0.20,0.20);
    lin += 2.0*bac*vec3(0.50,0.60,0.70);
    lin += 2.0*sss*vec3(0.20,0.20,0.20)*(0.2+0.8*dif*sha)*mate.w;
    lin += 2.0*spe*vec3(1.0)*sha*(0.3+0.7*fre);

    // surface-light interacion
    vec3 col = mate.xyz * lin;
    return col;
}

// fish type
#define FISH_WINE_SNAPPER   0
#define FISH_GREEN_SNAPPER  1
#define FISH_RED_SNAPPER    2
#define FISH_GOLD_SNAPPER   3
#define FISH_AMERICAN_SHAD  4

// (only for internal use)
vec3 fshPos = vec3(0.0);    // global coords
float fshTime = 0.0;        // fish motion time
vec2 fshCoord = vec2(0.0);  // fish length coord + fish width coord
vec3 fshSize = vec3(0.5);   // fish (width, height, length)
int fshType = 0;            // fish type ==> FishMaterial*

// Fish Definition
// 1. geometry (fishPos, fishTime, fishSize)
// 2. material (fishType)
#define FISH_AMOUNT         4
#define FISH_MAX_SIZE       10.0
vec3 fishPos[FISH_AMOUNT];  // fish position in world space
float fishTime[FISH_AMOUNT];// fish motion time
vec3 fishSize[FISH_AMOUNT]; // fish (width, height, length)
int fishType[FISH_AMOUNT];  // fish type ==> FishMaterial*

// fish material(color) <=== fish type
struct FishMaterial {
    float   shininess;      // 8.0;
    vec3    upperBodyCol;   // vec3(0.24,0.17,0.22);
    vec3    lowerBodyCol;   // vec3(1.0);
    float   upperSideLine;  // 0.4;
    float   lowerSideLine;  // 0.0;
    vec3    tailCombCol;    // vec3(2.0,1.0,0.5);
    vec3    tailPartCol;    // 0.9*vec3(2.6,1.5,1.0);
    vec3    stripeCol;      // vec3(0.5);
    float   bellyPattern;   // 50.0; (freq of belly pattern)
    vec3    topFinCol;      // vec3(0.8,0.2,0.2);
};
FishMaterial getFishMaterial( int type )
{
    FishMaterial fish;
    if( type == FISH_WINE_SNAPPER )
    {
        fish.shininess = 64.0;
    #ifdef ENABLE_FISH_COLORING
        fish.upperBodyCol = mix( vec3(0.3,0.17,0.17), vec3(0.17,0.17,0.17), 0.5+0.5*sin(time) );
    #else
        fish.upperBodyCol = vec3(0.3,0.17,0.17);
    #endif
        fish.lowerBodyCol = vec3(1.0);
        fish.upperSideLine = 0.4;
        fish.lowerSideLine = 0.0;
        fish.tailCombCol = 1.0*vec3(2.0,1.0,0.5);
        fish.tailPartCol = 0.9*vec3(2.6,1.5,1.0);
        fish.stripeCol = vec3(0.5);
        fish.bellyPattern = 50.0; // freq of belly pattern
        fish.topFinCol = vec3(0.8,0.2,0.2);
    }
    else if( type == FISH_GREEN_SNAPPER )
    {
        fish.shininess = 32.0;//8.0;
        fish.upperBodyCol = vec3(0.5,1.0,0.2);
        fish.lowerBodyCol = vec3(0.9,0.7,0.5);
        fish.upperSideLine = 0.2;
        fish.lowerSideLine = -0.6;
        fish.tailCombCol = vec3(1.0,1.0,0.0);
        fish.tailPartCol = vec3(1.0,1.0,0.0);
        fish.stripeCol = vec3(1.0,0.0,1.0);
        fish.bellyPattern = 75.0;//50.0;
        fish.topFinCol = vec3(1.5,0.0,0.0);
    }
    else if( type == FISH_RED_SNAPPER )
    {
        fish.shininess = 64.0;//8.0;
    #ifdef ENABLE_FISH_COLORING
        fish.upperBodyCol = mix( vec3(1.5,0.4,0.2), vec3(0.2,1.5,0.4), 0.5+0.5*sin(time) );
    #else
        fish.upperBodyCol = vec3(1.5,0.4,0.2);
    #endif
        fish.lowerBodyCol = vec3(0.9,0.7,0.5);
        fish.upperSideLine = 0.2;
        fish.lowerSideLine = -0.6;
        fish.tailCombCol = vec3(1.0,1.0,0.0);
        fish.tailPartCol = vec3(1.0,1.0,0.0);
        fish.stripeCol = vec3(1.5,0.1,0.1);
        fish.bellyPattern = 75.0;//50.0;
        fish.topFinCol = vec3(1.5,0.0,0.0);
    }
    else if( type == FISH_GOLD_SNAPPER ) // Goldfish
    {
        fish.shininess = 8.0;
        fish.upperBodyCol = vec3(0.9, 0.34, 0.07);
        fish.lowerBodyCol = vec3(0.85, 0.84, 0.81);
        fish.upperSideLine = 0.2;
        fish.lowerSideLine = -0.5;
        fish.tailCombCol = vec3(0.77, 0.48, 0.11);
        fish.tailPartCol = vec3(0.83, 0.38, 0.1);
        fish.stripeCol = vec3(0.88, 0.31, 0.15);
        fish.bellyPattern = 20.0;//50.0;
        fish.topFinCol = vec3(0.85, 0.13, 0.0);
    }
    else if( type == FISH_AMERICAN_SHAD )
    {
        fish.shininess = 8.0;
        fish.upperBodyCol = vec3(0.22, 0.35, 0.32);
        fish.lowerBodyCol = vec3(0.88, 0.91, 0.89);
        fish.upperSideLine = 0.4;
        fish.lowerSideLine = 0.1;
        fish.tailCombCol = vec3(0.65, 0.71, 0.6);
        fish.tailPartCol = vec3(0.58, 0.62, 0.32);
        fish.stripeCol = vec3(0.22, 0.35, 0.32);
        fish.bellyPattern = 10.0;//50.0;
        fish.topFinCol = vec3(0.42, 0.41, 0.33);
    }
    return fish;
}

vec2 vSegment( vec3 a, vec3 b, vec3 p )
// segment (has the length from a to b) 
// return.x = d * d (d = distance from p to segment)
// return.y = t (where a + t(b-a) = projection of p to segment) (0 < t < 1)
{
	vec3  pa = p - a;
	vec3  ba = b - a;
	float t = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	vec3  v = pa - ba*t;
	return vec2( dot(v,v), t );
}
float fishMap( in vec3 p )
// return dist from p to fish
// p ==> fishPos[i] ==> i ==> (fshPos, fshTime, fshSize, fshType, fshCoord)
{
    // vec3 spacing = vec3(3.0, 0.0, 7.0);
    // p = sRepeat( p, spacing );
    // p = sRotateY( p, PI*2.0 );

    //
    // find fishPos[i] which is the closest to p
    //
    float minDist = 1000.0;
    for(int i=0; i<FISH_AMOUNT; i++)
    {
        float dist = length( p - fishPos[i] );
        if( dist < minDist )
        {
            minDist = dist;
            fshPos = fishPos[i];
            fshTime = fishTime[i];
            fshSize = fishSize[i];
            fshType = fishType[i];
        }
    }

    // res.x = distance from p to fish
    // res.y = fish length coord
    // res.z = fish width coord
    vec3 res = vec3( 1000.0, 0.0, 0.0 );

	p -= fshPos;

    if( dot(p,p) > FISH_MAX_SIZE*FISH_MAX_SIZE ) return FISH_MAX_SIZE;

    float fishW = mapLinear( fshSize.x, 0.0, 1.0, 2.5, 0.8 ); // fish width:  2.5(thin) ~ 0.8(fat)
    float fishH = mapLinear( fshSize.y, 0.0, 1.0, 1.2, 0.5 ); // fish height: 1.2(short) ~ 0.5(tall)
    float fishL = mapLinear( fshSize.z, 0.0, 1.0, 2.0, 0.5 ); // fish length: 3.0(short) ~ 0.5(long)
    p *= vec3(fishW, fishH, fishL);

	vec3 q = p;
	
    vec3 a = vec3(0.0);
	a.x -= 0.25*sin(8.0 * 0.2 * fshTime);
	vec3 oa = a;

	float or = 0.0;
	float th = 0.0;
	float hm = 0.0;

	const int NUMI = 7;
	const float NUMF = 7.0;
	vec3 p1 = a; vec3 d1 = vec3(0.0);
	vec3 p2 = a; vec3 d2 = vec3(0.0);
	vec3 mp = a;
	for(int i=0; i<NUMI; i++)
	{	
		float ih = float(i)/NUMF;
		float ang = or + 1.0*(0.2+0.8*ih)*sin(3.0*ih - 2.0*fshTime);
		float ll = 0.26; if( i == (NUMI-1) ) ll = 0.4;
		vec3 b = a + ll*vec3(sin(ang), 0.0, cos(ang))*(16.0/NUMF);

		vec2 dis = vSegment( a, b, p );

		if( dis.x < res.x ) {
            res = vec3(dis.x, ih + dis.y / NUMF, 0.0);
            mp = a + (b - a)*dis.y;
        }

		if( i == 1 ) { p1 = a; d1 = b - a; }

		a = b;
	}
	float h = res.y;
	float ra = 0.04 + h*(1.0-h)*(1.0-h)*2.7;

	// tail
	p.y /= 1.0 + 14.0*(1.0-smoothstep( 0.0, 0.13, 1.0-h));
    p.z += 0.08*(1.0-saturate(abs(p.y)/0.075))*(1.0-smoothstep(0.0,0.1,1.0-h));
	res.x = 0.75 * (distance(p, mp) - ra);
	
	// mouth
	float d3 = 0.75*( length((p-oa)*vec3(0.5,2.0,1.0)) - 0.12 );
	res.x = max( -d3, res.x );

	// upper central fin
	float fh = smoothstep(0.15,0.2,h) - smoothstep(0.25,0.8,h);
	fh -= 0.2*pow(0.5+0.5*sin(210.0*h),0.2)*fh;
	d3 = length(p.xz-mp.xz) - 0.01;
    d3 = max( d3, p.y - (mp.y+ra+0.2*fh) );
	d3 = max( d3, -p.y - 0.0 );
	res.x = min( res.x, d3 );

	// fins
	d1.xz = normalize(d1.xz);

	float flap = 0.7 + 0.3*sin(2.0*8.0*0.2*fshTime);
    vec2 dd = normalize(d1.xz + sign((p-p1).x)*flap*d1.zx*vec2(-1.0,1.0));
	mat2 mm = mat2( dd.y, dd.x, -dd.x, dd.y );

	vec3 sq = p-p1;
	sq.xz = mm*sq.xz;
	sq.y += 0.2;
	sq.x += -0.15;
	float d = length( (sq-vec3(0.5,0.0,0.0))*vec3(1.0,2.0,1.0) ) - 0.3;
	d = 0.5*max( d, fBox( sq, vec3(1.0,1.0,0.01) ) );
    if( d < res.x ) res.z = smoothstep( 0.2, 0.7, sq.x );
	res.x = fsmin( d, res.x, 0.05 );

	sq = p-p1;
	sq.xz = mm*sq.xz;
	sq.y += 0.2;
	sq.x += 0.15;
	d = length( (sq-vec3(-0.5,0.0,0.0))*vec3(1.0,2.0,1.0) ) - 0.3;
	d = 0.5*max( d, fBox( sq, vec3(1.0,1.0,0.01) ) );
    if( d < res.x ) res.z = smoothstep( 0.2, 0.7, sq.x );
	res.x = fsmin( d, res.x, 0.05 );

    fshCoord = res.yz;

	return res.x;
}

vec3 fishColor( in vec3 rd, in vec3 pos, in vec3 nor, 
    in vec3 fshPos, in float fshTime, in vec2 fshCoord, in vec3 fshSize, in int fshType )
// fshPos, fshTime, fshCoord, fshSize ==> (geometry data)
// fshType ==> FishMaterial ==> (material(color) data)
{
    FishMaterial fish = getFishMaterial( fshType );

    // whole body
    vec4 mate = vec4(0.0);
    mate.xyz = fish.upperBodyCol; // material diffuse
    mate.w = fish.shininess;      // material shininess
    vec3 te = 0.8+2.2*texture( iChannel0, vec2(2.0*fshCoord.x,pos.y) ).xyz;
    mate.xyz *= te;
    
    // belly/backfin
    vec3 tailCol = mix( 1.0+0.5*sin(150.0*pos.y - sign(pos.y)*fshCoord.x*300.0), // <=== pattern of tail comb
                        1.0, smoothstep(0.0,0.1,1.0-fshCoord.x) ) * fish.tailPartCol;  // <=== tail part
    tailCol += (1.0-smoothstep(0.0,0.09,1.0-fshCoord.x)) * fish.tailCombCol;
    float iscola = smoothstep( 0.0, 0.2, 1.0-fshCoord.x );
    mate.xyz = mix( mix( vec3(te.x*0.5 + 1.5), tailCol, 1.0-iscola ) * 0.5, // <=== lower body
                    mate.xyz,                                               // <=== upper body
                    smoothstep(fish.lowerSideLine, fish.upperSideLine, nor.y) // upper if nor.y > upperSideLine, lower if nor.y < lowerSideLine
                );

    // stripes
    mate.xyz = mix( mate.xyz, 
                    (te.x+0.5)*1.0*fish.stripeCol, 
                    0.75*smoothstep( 0.5, 1.0, sin(1.0*te.x + fshCoord.x*100.0 + 13.0*nor.y) )*smoothstep(0.0,0.5,nor.y) );

    // escamas
    float ll = saturate( (fshCoord.x-0.2)/(0.8-0.2) );
    float ha = 1.0-4.0*ll*(1.0-ll);
    float pa = smoothstep( -1.0+2.0*ha, 1.0, sin(fish.bellyPattern*(pos.y-fshPos.y) ) ) * smoothstep( -1.0, 0.0, sin( 560.0*fshCoord.x ) );
    pa *= 1.0-smoothstep( 0.1, 0.2, nor.y );
    mate.xyz *= 0.5 + 0.5*fish.lowerBodyCol * (1.0-pa);

    // eye
    float r = length(vec2(5.0*fshCoord.x,pos.y)-vec2(0.5,0.13+fshPos.y) );
    r /= 1.2;
    mate.xyz = mix( mate.xyz, vec3(1.5)*saturate(1.0-r*4.0), 0.5*(1.0-smoothstep(0.08,0.09,r)) );
    mate.xyz *= smoothstep(0.03,0.05,r);
    mate.xyz += vec3(4.0)*(1.0-smoothstep(0.0,0.1,r))*
        pow( texture( iChannel1, 4.0*vec2(0.2*fshPos.z+4.0*fshCoord.x, pos.y) ).x, 2.0 );
    r = length(vec2(5.0*fshCoord.x,pos.y)-vec2(0.48,0.14) );
    mate.xyz = mix( mate.xyz, vec3(2.0), (1.0-smoothstep(0.0,0.02,r)) );
    
    // mouth
    vec3 oa = fshPos;
    oa.x -= 0.25*sin(8.0*0.2*fshTime);
    mate.xyz *= 0.1 + 0.9*step( 0.0, length( (pos - oa+vec3(0.0,0.0,-0.02))*vec3(1.5,2.0,1.0) ) - 0.14 );
    
    // top fin
    float fh = smoothstep(0.15,0.2,fshCoord.x) - smoothstep(0.25,0.8,fshCoord.x);
    float ra = 0.04 + fshCoord.x*(1.0-fshCoord.x)*(1.0-fshCoord.x)*2.7;
    float vv = saturate( (pos.y-fshPos.y-fshSize.y-ra-0.1) / 0.2 );
    vec3 fincol = mix( 1.0+0.5*sin(520.0*fshCoord.x), 1.0, vv ) * mix( fish.topFinCol, vec3(1.5,1.4,1.5), vv );
    mate.xyz = mix( mate.xyz, fincol, fh*smoothstep(0.0,0.05,pos.y-fshPos.y-fshSize.y-ra-0.1) );
    
    // side fins
    float isFin = fshCoord.y;
    fincol = 0.5*vec3(3.0,2.0,2.0) * mix(1.0+0.2*sin(150.0*pos.y),1.0,0.0);
    mate.xyz = mix( mate.xyz, fincol, isFin );
    mate.xyz *= 0.17;

    // lighting
    vec3 lig = SUN_LIGHT;
    float sky = saturate(nor.y);
    float bou = saturate(-nor.y);
    float dif = max(dot(nor,lig),0.0);
    float bac = max(0.3 + 0.7*dot(nor,-vec3(lig.x,0.0,lig.z)),0.0);
    float sha = 0.0; if( dif>0.001 ) sha = sceneShadow( pos+0.01*nor, lig, 0.0005, 32.0 );
    float fre = pow( saturate( 1.0 + dot(nor,rd) ), 5.0 );
    float spe = max( 0.0, pow( saturate( dot(lig, reflect(rd,nor)) ), mate.w ) ) * mate.w;
    float sss = pow( saturate( 1.0 + dot(nor,rd) ), 3.0 );
    
    // lights
    vec3 lin = vec3(0.0);
    float cc  = 0.55*texture( iChannel2, 1.8*0.02*pos.xz + 0.007*time*vec2( 1.0, 0.0) ).x;
          cc += 0.25*texture( iChannel2, 1.8*0.04*pos.xz + 0.011*time*vec2( 0.0, 1.0) ).x;
          cc += 0.10*texture( iChannel2, 1.8*0.08*pos.xz + 0.014*time*vec2(-1.0,-1.0) ).x;
    cc = 0.6*(1.0-smoothstep( 0.0, 0.025, abs(cc-0.4))) + 0.4*(1.0-smoothstep( 0.0, 0.150, abs(cc-0.4)));
    dif *= 1.0 + 2.0*cc;

    lin += 3.5*dif*vec3(1.00,1.00,1.00)*sha;
    lin += 3.0*sky*vec3(0.10,0.20,0.35);
    lin += 1.0*bou*vec3(0.20,0.20,0.20);
    lin += 2.0*bac*vec3(0.50,0.60,0.70);
    lin += 2.0*sss*vec3(0.20,0.20,0.20)*(0.2+0.8*dif*sha)*mate.w;
    lin += 2.0*spe*vec3(1.0)*sha*(0.3+0.7*fre);

    // material-light interacion
    vec3 col = mate.xyz * lin;
    return col;
}

//===========================================================

vec2 sceneMap( in vec3 p )
{
    vec2 res = vec2( skyMap(p), MATERIAL_SKY );
    res = dUnion( res, vec2( aquaTerrainMap(p), MATERIAL_TERRAIN ) );
    res = dUnion( res, vec2( fishMap(p), MATERIAL_FISH ) );
    return res;
}

void rayMinMax( in vec3 ro, in vec3 rd, out float tmin, out float tmax )
{
    tmin = 1.0;
    tmax = FAR;
    float max_height = SKY_HEIGHT;
    float t = (max_height - ro.y) / rd.y;
    if( t > 0.0 )
    {
        if( ro.y > max_height ) tmin = max( tmin, t );
        else                    tmax = min( tmax, t );
    }
    else
    {
        if( ro.y > max_height ) tmin = tmax = 1.0;
    }
}

vec3 render( in vec3 ro, in vec3 rd )
{
    float tmin, tmax;
	rayMinMax( ro, rd, tmin, tmax );

    vec2 tm = bisectMarching( ro, rd, tmin, tmax );
    if( tm.x > tmax ) tm.y = MATERIAL_SKY;

    vec3 col = vec3(0.0);

    if( tm.y == MATERIAL_SKY )
    {
        gl_FragDepth = 0.99;

        col = skyColor( SUN_LIGHT, rd, 0.0 );
        applyClouds( col, iChannel2, ro, rd );
	}
    else if( tm.y == MATERIAL_TERRAIN )
    {
        vec3 pos = ro + rd * tm.x;
        vec3 nor = sceneNormal( pos );
        gl_FragDepth = getFragDepth( pos );

        col = aquaTerrainColor( rd, pos, nor );
    }
    else if( tm.y == MATERIAL_FISH )
    {
        vec3 pos = ro + rd * tm.x;
        vec3 nor = sceneNormal( pos );
        gl_FragDepth = getFragDepth( pos );

        col = fishColor( rd, pos, nor, fshPos, fshTime, fshCoord, fshSize, fshType );
    }

    // fog under the water...
    col *= 0.65;
#ifdef ENABLE_AUTO_VIEW
    float density = 0.005;//0.01;
#else
    float density = 0.025;//0.2;
#endif
    //float t = max(0.0, tm.x-1.3);
    float t = max(0.0, tm.x);
    float h = 1.0-exp(-density*t);
    col = mix( col*(1.0-h), UNDER_WATER_COLOR, h );

    // water surface...
#if 1
    if( rd.y < 0.0 )
    {
        float waterHeight = 2.0;//2.0

        vec4 waterCol = waterColor( SUN_COLOR, SUN_LIGHT, WATER_COLOR, ro - vec3(0.0, waterHeight, 0.0), rd, vec4(0.0) );//cloudy = 0.0
        if( 0.0 < waterCol.w && waterCol.w < tm.x )
        {
            // foam on water...
            float t = (waterHeight-ro.y)/rd.y;
            vec2 uv = (ro + rd * t).xz;
            float sur = texture( iChannel3, 0.06*uv ).x;
            sur = smoothstep( 0.5, 1.0, sur )*0.5 + 0.5*sur*sur*smoothstep(0.2, 1.0, texture( iChannel2, 1.0*uv ).x);
            waterCol.rgb = mix( waterCol.rgb, vec3(2.5), 0.5*sur ); // foamCol = vec3(2.5)

            // sun specular...
            // float sunAmount = saturate( dot(SUN_LIGHT, reflect( rd, vec3(0.0,1.0,0.0) ) ) );
            // waterCol.rgb += 0.2*vec3(1.0,0.95,0.9)*pow(sunAmount,16.0);
            // waterCol.rgb += 0.5*vec3(1.0,0.95,0.9)*pow(sunAmount,96.0);

            col = mix(col, waterCol.rgb, saturate(1.1+rd.y));
        }
    }
#endif

    // post-processing
    col = pow( saturate(col), vec3(0.45) );
    col = mix( col, vec3(dot(col,vec3(0.333))), -0.5 );
	col = 0.5*col + 0.5*col*col*(3.0-2.0*col);
	col *= smoothstep( 0.0, 1.0, time );

	return col;
}

#ifdef ENABLE_AUTO_VIEW
mat3 cameraAutoView( in sampler2D tex, out vec3 ro, out vec3 rd )
{
    float curTime = 0.2*time;
	vec3 ta = fishPos[0] - vec3(2.0, 0.0, -2.0);
    ta.y += 2.0*sin(time);

    float r = 12.0;//8.0;
	ro = ta + vec3(r*sin(curTime), r, r*cos(curTime));

    float fl = 2.5;//2.0 1.2
    vec2 xy = (-resolution.xy + 2.0*gl_FragCoord.xy)/resolution.y;
    mat3 cam = cameraMatrix( ro, ta );
	rd = normalize( xy.x*cam[0] + xy.y*cam[1] + fl*cam[2] );
    return cam;
}
#endif

void main()
{
    // setup for fish animation...

#ifdef ENABLE_AUTO_VIEW
    float speedScale = 1.0;//0.7;
#else
    float speedScale = 0.0;
#endif

    fishTime[0] = time + 1.1;
    fishTime[1] = time + 2.2;
    fishTime[2] = time + 0.1;
    fishTime[3] = time + 1.5;

#ifdef ENABLE_FISH_COLORING
    fishPos[0] = vec3( 0.0, 0.0, -speedScale*fishTime[0] );
    fishPos[1] = vec3( 4.5, 1.0+0.5*sin(time), -speedScale*fishTime[1] );
    fishPos[2] = vec3( -5.0, 0.2, -speedScale*fishTime[2] );
    fishPos[3] = vec3( -9.5, 0.3+0.5*sin(time+0.5), -speedScale*fishTime[3] );

    fishSize[0] = vec3(0.5+0.3*sin(time), 0.5, 0.65+0.25*sin(time)); // fish (width, height, length)
    fishSize[1] = vec3(0.7, 0.0, 0.7);
    fishSize[2] = vec3(0.3, 0.5+0.3*sin(time), 0.8);
    fishSize[3] = vec3(0.5+0.3*sin(time), 0.5, 0.4+0.25*sin(time));
#else
    fishPos[0] = vec3( 0.0, 0.0, -speedScale*fishTime[0] );
    fishPos[1] = vec3( 4.5, 0.4, -speedScale*fishTime[1] );
    fishPos[2] = vec3( -5.0, 0.2, -speedScale*fishTime[2] );
    fishPos[3] = vec3( -9.5, 0.3, -speedScale*fishTime[3] );

    fishSize[0] = vec3(0.5, 0.5, 0.9); // fish (width, height, length)
    fishSize[1] = vec3(0.7, 0.0, 0.7);
    fishSize[2] = vec3(0.3, 0.7, 0.6);
    fishSize[3] = vec3(0.8, 0.5, 0.5);
#endif

    fishType[0] = FISH_WINE_SNAPPER;    // purple(wine)
    fishType[1] = FISH_GREEN_SNAPPER;   // green
    fishType[2] = FISH_RED_SNAPPER;     // red
    fishType[3] = FISH_GOLD_SNAPPER;    // gold

#ifdef ENABLE_AUTO_VIEW
    vec3 ro, rd;
    cameraAutoView( iChannel0, ro, rd );
#else
    @import ./ray;
#endif

    vec3 col = render( ro, rd );
    col = Vignetting( col, 0.5 );
    gl_FragColor = vec4( col, 1.0 );
}