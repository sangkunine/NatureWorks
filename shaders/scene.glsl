#ifndef RAYMARCH_SCENE
#define RAYMARCH_SCENE

//==============================================================================
// getMaterial() <=== terrain2FS.glsl
//==============================================================================

#define MATERIAL_SKY        0.0
#define MATERIAL_TERRAIN    1.0
#define MATERIAL_WATER      2.0
#define MATERIAL_TREES      3.0
#define MATERIAL_CLOUDS     4.0
#define MATERIAL_FISH       10.0
#define MATERIAL_TEXTURE0  20.0
#define MATERIAL_TEXTURE1  21.0
#define MATERIAL_TEXTURE2  22.0
#define MATERIAL_TEXTURE3  23.0

struct ShadeMaterial
{
    vec3 albedo; // base color
    float gloss; // gloss = 1 - roughness
    float metalness; // not used...
};

vec4 texCube( in sampler2D tex, in vec3 p, in vec3 n );

ShadeMaterial getMaterial( vec3 p, vec3 n, float m )
{
    ShadeMaterial mtl;
    if( m == MATERIAL_SKY )
    {
        mtl.albedo = vec3(0.2, 0.5, 0.85);
        mtl.gloss = 0.0;
    }
    else if( m == MATERIAL_WATER )
    {
        mtl.albedo = vec3(1.0);
        vec2 dp = max(abs(dFdx(p.xz)), abs(dFdy(p.xz)));
        float gloss = max(dp.x, dp.y);
        gloss = exp2( -gloss * 0.3 );
        mtl.gloss = gloss;
    }
    else if( m == MATERIAL_TEXTURE0 )
    {
        vec3 col = texCube( textureMaps[0], p, n ).xyz;
        mtl.albedo = col;
        mtl.gloss = 0.0;
    }
    else if( m == MATERIAL_TEXTURE1 )
    {
        vec3 col = texCube( textureMaps[1], p, n ).xyz;
        mtl.albedo = col;
        mtl.gloss = 0.0;
    }
    else if( m == MATERIAL_TEXTURE2 )
    {
        vec3 col = texCube( textureMaps[2], p, n ).xyz;
        mtl.albedo = col;
        mtl.gloss = 0.0;
    }
    else if( m == MATERIAL_TEXTURE3 )
    {
        vec3 col = texCube( textureMaps[3], p, n ).xyz;
        mtl.albedo = col;
        mtl.gloss = 0.0;
    }
    return mtl;
}

//==============================================================================
// Basic Functions: 
//      fresnelGloss(), getExtinction(), cheapCurvature()
//==============================================================================

vec3 fresnelGloss( vec3 n, vec3 v, vec3 F0, float gloss )
{
    float dotNV = max(0.0, dot(n, v));
    return F0 + (vec3(1.0) - F0) * pow(1.0 - dotNV, 5.0) * pow(gloss, 20.0);
}

vec3 getExtinction( float dist, float density, vec3 col )
// dist = traveled distance
// density = optical density (= absorption coeff) (eg: WATER_OPACITY*16.0 for water)
// col = extinctCol (= 1.0-vec3(0.5,0.4,0.1) for water)
{
    //return exp2( -dist * density * col ); // exp2(x) = 2^x
    return exp( -dist * density * col ); // exp2(x) = 2^x
}

float cheapCurvature( in vec3 p )
// here, used to darken the crevices(= crack)
{
    const float eps = 0.05, amp = 4.0, ampInit = 0.5;
    vec2 e = vec2(-1.0, 1.0)*eps; //0.05->3.5 - 0.04->5.5 - 0.03->10.->0.1->1.
    float t1 = sceneMap(p + e.yxx).x, t2 = sceneMap(p + e.xxy).x;
    float t3 = sceneMap(p + e.xyx).x, t4 = sceneMap(p + e.yyy).x;
    return saturate((t1 + t2 + t3 + t4 - 4.0*sceneMap(p).x)*amp + ampInit);
}

//==============================================================================
// Cook-Torrance PBR
//==============================================================================

float DistributionGGX(vec3 N, vec3 H, float roughness) {
    float a = roughness*roughness;
    float a2 = a*a;
    float dotNH = max(dot(N, H), 0.0);
    float denom = (dotNH*dotNH * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    return a2 / max(denom, 0.001); // prevent divide by zero for roughness=0.0 and dotNH=1.0
}
float GeometrySchlickGGX(float dotNV, float roughness) {
    float r = roughness + 1.0;
    float k = (r*r) / 8.0;
    return dotNV / (dotNV * (1.0 - k) + k);
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness) {
    float dotNV = max(dot(N, V), 0.0);
    float dotLN = max(dot(N, L), 0.0);
    float ggx2 = GeometrySchlickGGX(dotNV, roughness);
    float ggx1 = GeometrySchlickGGX(dotLN, roughness);
    return ggx1 * ggx2;
}
vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}
vec3 getCookShading( 
    in vec3 albedo, float metallic, float roughness, 
    in vec3 ld, in vec3 lc, 
    in vec3 n, in vec3 rd )
// ld = light direction
// lc = light color
// NOTE: we consider only the diffuse & specular of directional light
{
    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);
    vec3 radiance = lc;
    vec3 N = n;
    vec3 V = -rd;
    vec3 L = ld;
    vec3 H = normalize(V + L);
    float dotHV = max(dot(H, V), 0.0);
    float dotNV = max(dot(N, V), 0.0);
    float dotLN = max(dot(N, L), 0.0);
    float NDF = DistributionGGX( N, H, roughness );
    float G = GeometrySmith( N, V, L, roughness );
    vec3 F = fresnelSchlick( dotHV, F0 );
    vec3 specular = (NDF * G * F) / max(4.0 * dotNV * dotLN, 0.001);
    vec3 kD = (vec3(1.0) - F) * (1.0 - metallic);
    return (kD * albedo / PI + specular) * radiance * dotLN;
}
void getCookShading( 
    in vec3 albedo, float metallic, float roughness, 
    in vec3 ld, in vec3 lc, 
    in vec3 n, in vec3 rd, 
    out vec3 diffuse, out vec3 specular )
// ld = light direction
// lc = light color
// NOTE: we consider only the diffuse & specular of directional light
{
    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);
    vec3 radiance = lc;
    vec3 N = n;
    vec3 V = -rd;
    vec3 L = ld;
    vec3 H = normalize(V + L);
    float dotHV = max(dot(H, V), 0.0);
    float dotNV = max(dot(N, V), 0.0);
    float dotLN = max(dot(N, L), 0.0);
    float NDF = DistributionGGX( N, H, roughness );
    float G = GeometrySmith( N, V, L, roughness );
    vec3 F = fresnelSchlick( dotHV, F0 );
    vec3 spec = (NDF * G * F) / max(4.0 * dotNV * dotLN, 0.001);
    vec3 kD = (vec3(1.0) - F) * (1.0 - metallic);
    vec3 diff = kD * albedo / PI;
    vec3 lint = radiance * dotLN;
    diffuse = diff * lint;
    specular = spec * lint;
}
void getCookShading( vec3 albedo, float gloss, vec3 lc, vec3 ld, vec3 n, vec3 rd, out vec3 diffuse, out vec3 specular )
// NOTE*: fresnel not included (only called for internal use)
{
    vec3 v = -rd;
	vec3 h = normalize( v + ld );
	float dotLN = max(dot(ld, n), 0.0);
	float dotNV = max(dot(v, n), 0.0);
	float dotNH = max(dot(n, h), 0.0);
    vec3 lightWt = lc * dotLN;
    diffuse = albedo/PI * lightWt;

	float a = 1.0 - gloss;
	float a2 = a * a;
	float denom = dotNH * dotNH * (a2 - 1.0) + 1.0;
	float D = a2 / (PI * denom * denom);
	float k = a / 2.0;
    float vis1 = 1.0 / ((dotLN + 0.0001) * (1.0 - k) + k);
    float vis2 = 1.0 / ((dotNV + 0.0001) * (1.0 - k) + k);
    float G = vis1 * vis2;
	specular = (D * G) * lightWt;
}

//==============================================================================
// Blinn-Phong Shading
//==============================================================================

vec3 getPhongShading( in vec3 col, in float s, in vec3 ld, in vec3 lc, in vec3 p, in vec3 n, in vec3 rd )
// s = shininess (= 32, 64, 128, 256...)
// ld = light direction
// lc = light color
{
    #if 0
        float shadow = sceneShadow(p, ld);
        float ao = sceneAO(p, n);
        float atten = 1.0;
        float lfactor = atten * shadow * ao;
    #else
        float lfactor = 1.0;
    #endif

    vec3 v = -rd;
    vec3 h = normalize(ld + v);
    float amb = 0.0;//0.1
    float diff = max(dot(ld,n),0.0);
    float spec = pow(max(dot(n,h),0.0), s);

    return amb + (col*diff + spec) * lc * (lfactor);

    // < Blinn-Phong reflection >
    // col = ka*ia + { kd(L*N)*id + ks(R*V)^alpha*is } * (atten*shadow*ao)
    // col = ka*ia + { albedo*(L*N) + (N*H)^alpha } * lightCol * (atten*shadow*ao)
    // R*V (phong model) = N*H (blinn-phong model) (where H=L+V)
    // kd = albedo, ks = vec3(1.0)
    // atten = lightPower/(dist*dist)
    // alpha = shininess (eg: 16, 32, 64, 128)
    // i = incoming light, k = material
    // d = diffuse, s = specular, a = ambient
}

//==============================================================================
// Oren-Nayar Diffuse Lighting
//==============================================================================

float orenNayarDiffuse( in vec3 l, in vec3 n, in vec3 v, float r )
// return the diffuse light intensity (before multiplying albedo)
// l = lightDir (surf to light)
// n = surface normal at the sample point
// v = viewDir (surf to eye)
// r = surface roughness
{
    float r2 = r*r;
    float a = 1.0 - 0.5*(r2/(r2+0.57));
    float b = 0.45*(r2/(r2+0.09));

    float nl = dot(n, l);
    float nv = dot(n, v);

    float ga = dot(v-n*nv, n-n*nl);

	return max(0.0,nl) * (a + b*max(0.0,ga) * sqrt((1.0-nv*nv)*(1.0-nl*nl)) / max(nl, nv));
}

//==============================================================================
// Reflection (not including shadows & specular for high performance) (<=== riverFS.glsl)
//==============================================================================

vec3 skyColor( vec3 ld, vec3 rd, float horiz );

vec4 reflectRayColor( in vec3 lc, in vec3 ld, in vec3 rd, in vec3 p, in vec3 n, in float FAR )
// rd = 1st rays
// p = position at 1st hit
// n = normal at 1st hit
// return vec4(reflectCol, travelDist)
{
    vec3 rd2 = reflect( rd, n );

    vec2 hit = bisectMarching( p, rd2, 0.01, FAR );

    if( hit.x > FAR ) return vec4( skyColor( ld, rd2, 0.0 ), hit.x );

    // material
    vec3 p2 = p + rd2 * hit.x;
    vec3 n2 = sceneNormal( p2, EPS );
    ShadeMaterial mtl = getMaterial( p2, n2, hit.y );

    // lighting
    vec3 diffuse, specular;
    getCookShading( mtl.albedo, mtl.gloss, lc, ld, n2, rd2, diffuse, specular );
    return vec4( diffuse, hit.x ); // we ignore shadows, specular...
}

//==============================================================================
// Refraction (or transmittance) with inscatter & extinction (<=== riverFS.glsl)
//==============================================================================

vec4 refractRayColor( in vec3 lc, in vec3 ld, in vec3 rd, in vec3 p, in vec3 n, in float eta, in float density, in float FAR )
// rd = 1st rays
// p = position at 1st hit
// n = normal at 1st hit
// eta = airIOR / waterIOR (= 1.0 / 1.3333) when light travel from air to water
// density = optical density (eg: WATER_OPACITY*6.0 for water)
// return vec4(refractCol, travelDist)
{
    vec3 rd2 = refract( rd, n, eta );

    vec2 hit = bisectMarching( p, rd2, 0.01, FAR );

    if( hit.x > FAR ) return vec4( skyColor( ld, rd2, 0.0 ), hit.x );

    // material
    vec3 p2 = p + rd2 * hit.x;
    vec3 n2 = sceneNormal( p2, EPS );
    ShadeMaterial mtl = getMaterial( p2, n2, hit.y );

    // lighting
    vec3 diffuse, specular;
    getCookShading( mtl.albedo, mtl.gloss, lc, ld, n2, rd2, diffuse, specular );
    // here, we ignore shadows, specular...

    // inscatter & extinction
    float travelDist = hit.x;
    float sunAmount = dot(ld, rd);
    vec3 inscatter = diffuse * (1.0 - exp( -travelDist * 0.1 )) * (1.0 + sunAmount);
    #if 0
        vec3 extinction = getExtinction( travelDist, density, diffuse ); // study required...
    #else
        vec3 extinction = getExtinction( travelDist, density, 1.0-vec3(0.5,0.4,0.1) );
    #endif

    return vec4( (diffuse + inscatter) * extinction, hit.x );
}

//==============================================================================
// flashColor() <=== weatherFS.glsl
//==============================================================================

#ifdef USE_FLASH_1
vec3 flashColor( in float cloudy, in float time )
{
    float lightning = 0.0;
    if( cloudy > 0.77 ) {
        float f = mod(time + 1.5, 2.5);
        if( f < 0.8 ) {
            f = smoothstep(0.8, 0.0, f) * 1.5;
            lightning = mod(-time*(1.5 - hash11(time*0.3)*0.002), 1.0) * f;
        }
    }
    return saturate(vec3(1.0, 1.0, 1.2) * lightning);
}
#endif

//==============================================================================
// applyClouds() <=== cavesFS.glsl + terrainFS.glsl + terrain2FS.glsl
//==============================================================================

void applyClouds( inout vec3 col, in vec3 ro, in vec3 rd )
{
    vec3 cloud = vec3(1.0, 0.95, 1.0);
    float amount = 300.0; // 100 ~ 500
    vec2 p = ro.xz + (rd.xz/rd.y) * (amount - ro.y);
    col = mix( col, cloud, 0.5*smoothstep(0.5, 0.8, fbm_4(0.005*p)) );
}

void applyClouds( inout vec3 col, in sampler2D tex, in vec3 ro, in vec3 rd )
{
    vec3 cloud = vec3(1.0, 0.95, 1.0);
    float amount = 300.0; // 100 ~ 500
    vec2 p = ro.xz + (rd.xz/rd.y) * (amount - ro.y);
    col = mix( col, cloud, 0.5*smoothstep(0.5, 0.8, fbm_4(tex, 0.005*p)) );
}

//==============================================================================
// applySnow()
//==============================================================================

void applySnow( inout vec3 col, in float hmin, in float hmax, in vec3 p, in vec3 n )
{
    // cf: soilColor = 0.1*vec3(0.62, 0.65, 0.7)
    // amount (based on height): smoothstep( h1, h2, y + random )
    // dense (based on normal): denser as n.y point toward sky
    vec3 snow = 0.28*vec3(0.62, 0.65, 0.7);
    //float amount = smoothstep( 55.0, 80.0, p.y + 25.0*fbm12(0.01*p.xz) );
    float amount = smoothstep( hmin, hmax, p.y + (hmax-hmin)*fbm12(0.01*p.xz) );
    float dense = smoothstep( 1.0-0.5*amount, 1.0-0.1*amount, n.y );
    float strength = amount*dense;
    col = mix( col, snow, smoothstep( 0.1, 0.9, strength ) );
}
void applySnow( inout vec3 col, in sampler2D tex, in float hmin, in float hmax, in vec3 p, in vec3 n )
// (hmin, hmax) = (55.0, 80.0)
{
    // cf: soilColor = 0.1*vec3(0.62, 0.65, 0.7)
    // amount (based on height): smoothstep( h1, h2, y + random )
    // dense (based on normal): denser as n.y point toward sky
    vec3 snow = 0.28*vec3(0.62, 0.65, 0.7);
    //float amount = smoothstep( 55.0, 80.0, p.y + 25.0*fbm(tex, 0.01*p.xz) );
    float amount = smoothstep( hmin, hmax, p.y + (hmax-hmin)*fbm(tex, 0.01*p.xz) );
    float dense = smoothstep( 1.0-0.5*amount, 1.0-0.1*amount, n.y );
    float strength = amount*dense;
    col = mix( col, snow, smoothstep( 0.1, 0.9, strength ) );
}
float snowFlake(vec2 uv, vec2 center, float radius)
{
    return 1.0 - sqrt(smoothstep(0.0, radius, length(uv - center)));
}
vec3 snowColor( float flakeSize, float windSpeed )
// snow falls from the sky (eg: col += snowColor(0.2, 0.5))
// flakeSize = size of snow flake (0.0 ~ 1.0, default=0.2)
// windSpeed = speed of wind (0.0 ~ 1.0, default=0.5)
{
	const float NUM_SHEETS = 10.0;
	const float NUM_FLAKES = 400.0;

	flakeSize *= 0.002;
	windSpeed *= 2.0;
	float windTime = windSpeed*time;
    vec2 uv = gl_FragCoord.xy / resolution.x;
	vec3 col = vec3(0.0);

    for(float i=1.0; i<=NUM_SHEETS; i+=1.0)
	{
        for(float j=1.0; j<=NUM_FLAKES; j+=1.0)
		{
            if( j > NUM_FLAKES / i ) break;
			float size = flakeSize*i * (1.0 + rand(j)/2.0);
            float speed = 0.75*size + rand(i)/1.5;

            vec2 center = vec2(0.0);
            center.x = -0.3 + rand(j*i) * 1.4 + 0.1*cos(windTime + sin(j*i));
            center.y = fract(sin(j) - speed * windTime) / 1.3;

            col += vec3( (1.0 - i/NUM_SHEETS) * snowFlake(uv, center, size) );
        }
    }
	return col;
}

//==============================================================================
// sunScatter()
//==============================================================================

vec3 sunScatter( vec3 ld, vec3 rd )
// sun glare...
{
    float sunAmount = max( dot(ld, rd), 0.0 );
    return 0.3*vec3(1.0,0.7,0.3) * pow( sunAmount, 8.0 );
    //return 0.25*vec3(1.0,0.4,0.2)*pow( sunAmount, 4.0 );
}

//==============================================================================
// skyColor()
//==============================================================================

vec3 skyColor( vec3 rd )
// sky only
{
    rd.y = max( rd.y, 0.0 );
    return vec3( pow(1.0-rd.y, 2.0), 1.0-rd.y, 0.6+(1.0-rd.y)*0.4 );
}
vec3 skyColor( vec3 ld, vec3 rd )
// ld = light direction (eg: SUN_LIGHT)
{
	float sunAmount = max( dot(ld, rd), 0.0 );
#if 0
	vec3 hor = mix( 1.2*vec3(0.7,1.0,1.0), vec3(1.5,0.5,0.05), 0.25+0.75*sunAmount );
	vec3 col = mix( vec3(0.2,0.6,0.9), hor, exp(-(4.0+2.0*(1.0-sunAmount))*max(0.0,rd.y-0.1)) );
#else
	//vec3 hor = vec3(1.2);
	vec3 hor = 2.5*vec3(0.48, 0.49, 0.53);//2.5*fogColor
	vec3 col = mix( vec3(0.2,0.3,0.9), hor, exp(-(4.0+2.0*(1.0-sunAmount))*(rd.y-0.1)) );
#endif
    col *= 0.5;
	col += 0.35*vec3(1.0,0.8,0.7)*pow(sunAmount,512.0);
	col += 0.25*vec3(1.0,0.4,0.6)*pow(sunAmount,32.0);
	col += 0.25*vec3(1.0,0.2,0.4)*pow(sunAmount,4.0);
	return col;
}
vec3 skyColor( vec3 lc, vec3 ld, vec3 rd )
// sky + sun(glare,sizable)
{
    vec3 sky = skyColor( rd );
    float sunAmount = max( dot(rd, ld), 0.0 );
    float sunGlare = 0.32; //[0.2, 0.4]
    float sunSize = 1.0;   //[0.3, 2.0]
	sky += lc * pow(sunAmount, 6.5) * sunGlare;
    sky += lc * min( pow( sunAmount, 1150.0 ), 0.7 ) * sunSize;
	return sky;
}
vec3 skyColorGradient( vec3 lc, vec3 ld, vec3 rd )
// sky + sun(glare)
{
    // sky
	float v = pow(1.0 - max(rd.y, 0.0), 5.0)*0.5;
	vec3 sky = v * lc*vec3(0.4) + vec3(0.18,0.22,0.4);
	// sun: wide glare effect...
    float sunAmount = max(dot(rd, ld), 0.0);
	sky += lc * min(pow(sunAmount, 60.5)*0.32, 0.3);
	sky += lc * min(pow(sunAmount, 1150.0), 0.3)*0.65;
	return sky;
}
vec3 skyColorGlare( vec3 lc, vec3 ld, vec3 rd )
{
    float sunAmount = max( dot(rd, ld), 0.0 );
	float v = pow(1.0 - max(rd.y,0.0), 6.0);
	vec3  sky = mix(vec3(0.1, 0.2, 0.3), vec3(0.32, 0.32, 0.32), v);
	sky += lc * sunAmount * sunAmount * 0.25;
	sky += lc * min(pow(sunAmount, 800.0)*1.5, 0.3);
	return saturate(sky);
}
vec3 skyColor( vec3 ld, vec3 rd, float horiz )
// sky + sun(layered)
// ld = direction from surface to sun
// horiz = 0.0(not added) or 1.0(add horiz)
{
    float sunAmount = max( dot(ld, rd), 0.0 );

    // sky
    vec3 skyTop = vec3(0.2, 0.4, 0.85); skyTop = 1.2*skyTop - rd.y*rd.y*0.5;
    vec3 skyBottom = vec3(0.6, 0.64, 0.72);
    vec3 sky = mix( skyTop, skyBottom, pow(1.0-max(rd.y,0.0), 4.0) );
    vec3 sunRedish = vec3(0.4, 0.375, 0.35);
    sky = mix( sky, sunRedish, sunAmount*0.75 );// add redish around the sun

    // sun with three layers (much better)
    vec3 sun = 0.25*vec3(1.0,0.7,0.4) * pow(sunAmount, 5.0);
    sun += 0.25*vec3(1.0,0.8,0.6) * pow(sunAmount, 64.0);
    sun += 0.2*vec3(1.0,0.9,0.7) * pow(sunAmount, 512.0);

    // sky with sun
    sky += sun;

    // sky with horizon
    if( horiz == 1.0 )
    {
        vec3 horizCol = 0.68*vec3(0.4, 0.65, 1.0);
        sky = mix( sky, horizCol, pow( 1.0-max(rd.y,0.0), 16.0 ) );
    }

    return sky;
}

//==============================================================================
// skyCloudsColor() <=== forestFS.glsl
//==============================================================================

vec3 skyCloudsColor( in sampler2D tex, in vec3 lc, in vec3 ld, in vec3 ro, in vec3 rd )
// tex = noise texture for clouds (eg: grayNoise256.png)
{
    // sky (background)
    vec3 col = 0.9*vec3(0.4,0.65,1.0) - rd.y*vec3(0.4,0.36,0.4);

    // clouds
    float t = (500.0 - ro.y) / rd.y;
    if( t > 0.0 )
    {
        vec2 uv = (ro + t*rd).xz;
        float cl = -1.0 + 2.0*fbm_9( tex, uv*0.002 );
        float dl = smoothstep(-0.2,0.6,cl);
        col = mix( col, vec3(1.0), 0.4*dl );
    }

	// sun glare
    float sunAmount = max( dot(ld, rd), 0.0 );
    col += 0.6*lc*pow( sunAmount, 32.0 ); // eg: lc = vec3(1.0,0.6,0.3)

	return col;
}

//==============================================================================
// applyFog() <=== pbrFS.glsl + riverFS.glsl + terrainFS.glsl + terrain2FS.glsl + forestFS.glsl
//==============================================================================

void applyFog( inout vec3 col, in vec3 fog, in float density, in float dist )
{
    float amount = 1.0 - exp(-pow(dist * density, 1.5));
    col = mix( col, fog, amount );
}
void applyFog( inout vec3 col, float dist )
{
    float density = 0.0005; // 0.0002 ~ 0.001
    vec3 fogCol = 0.65*vec3(0.4, 0.65, 1.0);
    applyFog( col, fogCol, density, dist );
}
void applyFog( inout vec3 col, in vec3 lc, in vec3 ld, in vec3 rd, in float density, in float dist )
// fog density ===> constant
{
    float fogAmount = 1.0 - exp(-dist * density);
    vec3 fogCol = skyColorGradient( lc, ld, rd );
    #ifdef USE_LIGHT_FLARE
        float dotLV = saturate(dot(ld, -rd));
        fogCol += lc * pow(dotLV, 10.0);
    #endif
    col = mix(col, fogCol, fogAmount);
}
void applyFog( inout vec3 col, in vec3 lc, in vec3 ld, in vec3 ro, in vec3 rd, in float dist )
// fog density ==> changed along y (cf: d(y) = a * e ^ (-b*y))
{
#if 1
    float c = 0.05;  // fog factor
    float b = 0.25;  // fog falloff
#else
    float c = 0.15;  // fog factor
    float b = 0.025; // fog falloff
#endif

    float fogAmount = c * exp(-ro.y * b) * (1.0 - exp(-dist * rd.y * b)) / rd.y;
    vec3 fogCol = skyColorGradient( lc, ld, rd );
    #ifdef USE_LIGHT_FLARE
        float dotLV = saturate(dot(ld, -rd));
        fogCol += lc * pow(dotLV, 10.0);
    #endif
    col = mix(col, fogCol, fogAmount);
}

//==============================================================================
// applyLensFlares() <=== weatherFS.glsl
//==============================================================================

#ifdef USE_LENSFLARES_1

void applyLensFlares( inout vec3 col, mat3 camMat, vec2 xy, vec3 lc, vec3 ld, float cloudy )
// camMat[0] = camera rightVec
// camMat[1] = camera upVec
// camMat[2] = camera forwardVec
// xy = [-1.7,1.7] x [-1,1] (1.77 = aspect)
// lc = light(sun) color, ld = light direction
{
	vec3 cu = camMat[0];
	vec3 cv = camMat[1];
	vec3 cw = camMat[2];
	float bri = dot(cw, ld) * 2.7 * clamp(-(0.7*cloudy-0.45) + 0.2, 0.0, 0.2);
	if( bri > 0.0 )
	{
		vec2 sunPos = vec2( dot( ld, cu ), dot( ld, cv ) );
		vec2 uvT = xy - sunPos;
		uvT = uvT * length( uvT );
		bri = pow( bri, 6.0 )*0.6;

		float glare1 = max(1.2 - length(uvT + sunPos*2.0)*2.0, 0.0);
		float glare2 = max(1.2 - length(uvT + sunPos*0.5)*4.0, 0.0);
		uvT = mix (uvT, xy, -2.3);
		float glare3 = max(1.2 - length(uvT + sunPos*5.0)*1.2, 0.0);

		col += bri * lc * vec3(1.0, 0.5, 0.2)  * pow(glare1, 10.0)*25.0;
		col += bri * vec3(0.8, 0.8, 1.0) * pow(glare2, 8.0)*9.0;
		col += bri * lc * pow(glare3, 4.0)*10.0;
	}
}
void applyLensFlares( inout vec3 col, in vec3 camPos, mat4 worldMat, mat4 iprojMat, vec2 xy, vec3 lc, vec3 ld, float cloudy )
// camPos = camera position (= ro)
// worldMat = view space to world space
// iprojMat = proj space to view space
// xy = [-1.7,1.7] x [-1,1] (1.77 = aspect)
// lc = light(sun) color, ld = light direction
// cloudy = cloud amount (0 ~ 1)
{
    vec3 camTar = (worldMat * iprojMat * vec4(0.,0.,0.,1.)).xyz;
	mat3 camMat = cameraMatrix(camPos, camTar);
    applyLensFlares( col, camMat, xy, lc, ld, cloudy );
}
#endif

//==============================================================================
// cloudsColor() (volumetric) <=== weatherFS.glsl
//==============================================================================

#ifdef USE_CLOUDS_1

const float CLOUD_LOWER = 2800.0;
const float CLOUD_UPPER = 3800.0;

float cloudsMap( vec3 p, float cloudy )
{
	float h = -(fbm13(p*0.0005) - (0.7*cloudy-0.45) - 0.6);
	return h;
}
float cloudsMapBounded( vec3 p, float cloudy )
{
	float h = cloudsMap( p, cloudy );
    h *= smoothstep(CLOUD_UPPER + 100.0, CLOUD_UPPER, p.y);
	return h;
}
float cloudsInScattering( vec3 p, vec3 ld, float cloudy )
// return intensity of light incident on p(inside the cloud) along the ray(ld) from sun
// where ld = ray direction from p to sun
{
	float cloudThick = CLOUD_UPPER - CLOUD_LOWER;
	cloudThick *= 0.2;
	//float cloudThick = 1.0;
    float l = cloudsMapBounded( p, cloudy ) - cloudsMapBounded( p + ld * cloudThick, cloudy );
	return clamp( -l*2.0, 0.05, 1.0 );
}
vec3 cloudsColor( in vec3 lc, in vec3 ld, in vec3 ro, in vec3 rd, in vec4 cloudy, out float density )
// lc = light color, ld = light direction
// cloudy.x = cloud amount (0.0 ~ 1.0)
// cloudy.yzw = cloud flash color
// density = clouds density accumulated along the ray(rd)
{
	float tmin = ((CLOUD_LOWER-ro.y) / rd.y);
	float tmax = ((CLOUD_UPPER-ro.y) / rd.y);
	vec3 p = ro + rd*tmin;
	#if 1
		tmin += hash12(p.xz*19.19)*350.0;
	#endif

	// trace the inside of clouds
	const int iter = 20;//55
	vec3 dp = rd * (tmax-tmin) / float(iter);
	vec2 shadeSum = vec2(0.0, 0.0);
	for(int i = 0; i < iter; i++)
	{
		if( shadeSum.y >= 1.0 ) break;
        vec2 shade;
		shade.x = cloudsInScattering(p, ld, cloudy.x); // shade.x = light intensity
		shade.y = max(cloudsMap(p, cloudy.x), 0.0);    // shade.y = cloud density
		//shade.x *= shade.y;
		shadeSum += shade * (1.0 - shadeSum.y);        // accumulation
		p += dp;
	}
	vec3 clouds = mix(vec3(pow(shadeSum.x, 0.6)), lc, (1.0-shadeSum.y)*0.4);

	// add flash (= cloudy.yzw)
    clouds += cloudy.yzw * (shadeSum.y + shadeSum.x + 0.2) * 0.5;
	density = shadeSum.y;
	return clouds;
}

//==============================================================================
// skyCloudsColor() <=== weatherFS.glsl
//==============================================================================

vec3 skyCloudsColor( in vec3 lc, in vec3 ld, in vec3 ro, in vec3 rd, in vec4 cloudy )
// skyColor() + cloudsColor()
// ld = light direction toward the sun
// cloudy.x = cloud amount (0.0 ~ 1.0)
// cloudy.yzw = cloud flash color
{
	float density;
	vec3 sky = skyColor( ld, rd, 1.0 );
	vec3 clouds = cloudsColor( lc, ld, ro, rd, cloudy, density );
	sky = mix( sky, min(clouds, 1.0), density );
	return saturate( sky );
}
#endif

//==============================================================================
// cloudsColor() (volumetric) <=== forestFS.glsl
//==============================================================================

#ifdef USE_CLOUDS_2

const float CLOUD_LOWER_HEIGHT = 50.0;
const float CLOUD_UPPER_HEIGHT = 500.0;
const float CLOUD_AMOUNT = 0.0001; // 0.0 ~ 1.0
const vec3 FOG_COLOR = vec3(0.4, 0.6, 1.15);

vec4 cloudsMap( in vec3 pos )
// return (x,yzw) = (density, gradient)
{
    vec4 n = fbmd_7( pos*0.003*vec3(0.6,1.0,0.6) - vec3(0.1,1.9,2.8) );
    n.x = -1.0 + 2.0*n.x;
    n.yzw = 2.0 * n.yzw;
    float h0 = CLOUD_LOWER_HEIGHT;
    float h1 = CLOUD_UPPER_HEIGHT;
    float hm = mix( h0, h1, 0.1 );
    vec2 h =  smoothstepd( h0, hm, pos.y ) -  smoothstepd( hm, h1, pos.y );
    h.x = 2.0*n.x + h.x + mapLinear( CLOUD_AMOUNT, 0.0, 1.0, -1.5, 0.0 );
    return vec4( h.x, 2.0*n.yzw*vec3(0.6,1.0,0.6)*0.003 + vec3(0.0,h.y,0.0) );
}
vec4 cloudsColor( in vec3 ld, in vec3 ro, in vec3 rd, in float tmin, in float tmax )
// ld = light direction (eg: SUN_LIGHT)
{
    vec4 sum = vec4(0.0);

    // bounding volume!!
    float tl = (CLOUD_LOWER_HEIGHT - ro.y) / rd.y;
    float th = (CLOUD_UPPER_HEIGHT - ro.y) / rd.y;

    tl = max( tl, 0.0 );
    th = max( th, 0.0 );
    tmin = max( tmin, min( tl, th ) );
    tmax = min( tmax, max( tl, th ) );

    float t = tmin;
    float thickness = 0.0;
    float delta = (tmax - tmin)/128.0;

    for(int i=0; i<128; i++)
    {
        vec3  pos = ro + t*rd;
        vec4  denGrad = cloudsMap( pos ); 
        float den = denGrad.x;
        float dt = max(delta, 0.011*t);//0.1

        if( den > 0.001 )
        {
        #if 1
            float sha = 1.0; // low quality
        #else
            float sha = clamp( 1.0 - max(0.0, cloudsMap( pos + ld*5.0 ).x), 0.0, 1.0 );
        #endif

            // lighting
            vec3 n = -normalize( denGrad.yzw );
            float dif = saturate( dot(n, ld) )*sha; 
            float fre = saturate( 1.0 + dot(n,rd) )*sha;
            vec3 lin  = vec3(0.70,0.80,1.00)*0.9*(0.6+0.4*n.y);
            lin += vec3(0.20,0.25,0.20)*0.7*(0.5-0.5*n.y);
            lin += vec3(1.00,0.70,0.40)*4.5*dif*(1.0-den);        
            lin += vec3(0.80,0.70,0.50)*1.3*pow(fre,32.0)*(1.0-den);

            // color
            vec3 col = vec3(0.8,0.77,0.72) * saturate(1.0-4.0*den);
            col *= lin;

            // fog added
            applyFog( col, FOG_COLOR, 0.001, t );

            // front to back blending
            float alpha = saturate( den*0.25*min(dt, tmax-t-dt) );
            col.rgb *= alpha;
            sum = sum + vec4(col, alpha)*(1.0 - sum.a);

            thickness += dt * den;
        }
        else
        {
            dt *= 1.0 + 4.0*abs(den);
        }
        t += dt;
        if( sum.a > 0.995 || t > tmax ) break;
    }
    
    if( thickness > 0.0 )
    {
        float sunAmount = saturate( dot(ld, rd) );
		sum.xyz += vec3(1.0,0.6,0.4)*0.2 * pow(sunAmount,32.0) * exp(-0.3*thickness) * saturate(thickness*4.0);
    }

    return saturate( sum );
}
#endif

//==============================================================================
// applyRain() <=== weatherFS.glsl
//==============================================================================

#ifdef USE_RAIN_1
void applyRain( inout vec3 col, in sampler2D tex, in vec2 xy, in vec4 cloudy, in float time )
// xy = [-1.7,1.7] x [-1,1] (1.77 = aspect)
// cloudy.x = cloud amount (0.0 ~ 1.0)
// cloudy.yzw = cloud flash color
{
    vec2 uv = gl_FragCoord.xy / resolution.xy;
	vec2 st = xy * vec2(0.5+(uv.y+1.0)*0.3, 0.02) + vec2(time*0.5+uv.y*0.2, time*0.2);
 	float f = texture(tex, st, -100.0).y * texture(tex, st*0.773, -100.0).x * 1.55;
	float rain = saturate((0.7*cloudy.x-0.45) - 0.15);
	f = clamp( pow(abs(f), 15.0)*5.0*(rain*rain*125.0), 0.0, (uv.y+0.1)*0.6 );
	col = mix( col, vec3(0.15) + cloudy.yzw, f );
	col = saturate(col);
}
#endif

//==============================================================================
// waterColor() for simple/fast waters (under sky & clouds) <=== weatherFS.glsl + terrain2FS.glsl
//==============================================================================

#if defined(USE_WATER_1) || defined(USE_WATER_2)
const float WATER_HEIGHT = 0.0;
const float WATER_SPEED = 0.5;      // 0.0 ~ 1.0
const float WATER_CHOPPINESS = 0.5; // 0.0 ~ 1.0 (not used...)
const float WATER_OPACITY = 0.5;    // 0.0 ~ 1.0
#endif

#ifdef USE_WATER_1
vec4 waterColor( in vec3 lc, in vec3 ld, in vec3 wc, in vec3 ro,in vec3 rd, in vec4 cloudy )
// lc = light(sun) color, ld = light direction
// wc = water color (eg: vec3(0.3, 0.4, 0.45))
// cloudy.x = cloud amount (0.0 ~ 1.0)
// cloudy.yzw = cloud flash color
// reflected: skyColor() + cloudsColor(cloudy)
// if 0.0 < vec4.w < sceneDist, then ray hits the water before hitting the terrain scene (see terrain2FS.glsl)
{
	// water position (where ray intersects with water)
	float t = (WATER_HEIGHT-ro.y)/rd.y;
    if( t < 0.0 ) return vec4(0.0,0.0,0.0,t);
	vec3 p = ro + rd * t; // p = waterPos

	// water normal
	float distort = WATER_SPEED * time;
	float tx = cos(p.x*.052)*4.5;
	float tz = sin(p.z*.072)*4.5;
	vec2 co = noise22( vec2(p.x*4.7 + 1.3 + tz, p.z*4.69 + distort*35.0 - tx) );
	co += noise22( vec2(p.z*8.6 + distort*13.0 - tx, p.x*8.712 + tz) )*0.4;
	vec3 n = normalize( vec3(co.x, 20.0, co.y) );

	// reflected: sky + clouds
	vec3 re = reflect(rd, n);
    #ifdef USE_CLOUDS_1
        vec3 sky = skyCloudsColor(lc, ld, p, re, cloudy);
    #else
        float density;
        vec3 sky = skyColor( ld, re, 1.0 );
        applyClouds( sky, ro, rd );
        sky = saturate( sky );
    #endif

	// lighting...
	#if 1
		float fresnel = max(dot(n, -rd), 0.0);
		fresnel = pow(fresnel, 0.3) * 1.1;
		vec3 water = mix( wc * max(dot(ld, n), 0.0), sky*0.6, fresnel );
	#else
		// mix sea color (diffuse) with sky color (specular)
		float FAR = CLOUD_UPPER;
		float atten = smoothstep( FAR, FAR*0.7, t );
		float F0 = 0.0204; // water
		vec3 V = -rd;
		vec3 H = normalize(ld + V);
		float F = mix(F0, 1.0, pow(1.0 - max(dot(H, V),0.0), 5.0));// schlick approximation
		vec3 water = mix( atten*wc*max(dot(ld, n), 0.0), sky*0.6, F );
	#endif

	// optional: add the reflected sun(= lc)...
	#if 1
		float glit = max(dot(re, ld), 0.0);
		water += lc * pow(glit, 220.0) * max(-(0.7*cloudy.x-0.45)*100.0, 0.0);
	#endif

	// optional: add the water glint...
	#if 1
		tx = p.y - ro.y;
		water += vec3(0.1)*saturate( 1.0 - pow(tx + 0.5, 3.0) * texture(textureMaps[0], p.xz*0.1, -2.0).x );
	#endif

	return vec4(water, t);
}
#endif

//==============================================================================
// waterColor() for river waters (under sky + above terrain) <=== riverFS.glsl
//==============================================================================

#ifndef NEW_RIVER_DEPTH
const float RIVER_WAVE_LENGTH = 16.0;
const float RIVER_WAVE_HEIGHT = 1.5;
float riverCurve( float x ) {
    // definition of meandering river
    float L = RIVER_WAVE_LENGTH;
    float H = RIVER_WAVE_HEIGHT;
    return H * sin( PI2/L * x );
}
float riverCurved( float x ) {
    // derivative of riverCurve() above
    float L = RIVER_WAVE_LENGTH;
    float H = RIVER_WAVE_HEIGHT;
    float W = PI2/L;
    return H * W * cos(W * x);
}
float riverBaseDepth( vec3 p ) {
    float depth = 0.3 + (0.5 + 0.5*sin(p.x*0.001 + 3.0)) * 0.4;
    float width = 2.0 + cos( p.x * 0.1 ) * 1.0;
    float amount = smoothstep( width, width * 0.5, abs(p.z - riverCurve(p.x)) );
    return amount * depth;
}
float riverDepth( vec3 p ) {
    float depth = 0.0;
    depth += riverBaseDepth(p);
    depth += 0.45*riverBaseDepth(p*2.0);
    return depth;
}
#endif

#ifdef USE_WATER_2
vec3 waterNoise( vec2 waterPos, vec2 flowOffset, float foamScale, float gradAscent )
// foamScale = scaling factor of noise amplitude
// gradAscent = scaling factorof noise derivative
{
    waterPos *= (WATER_CHOPPINESS*2.0);
	vec3 f = vec3(0.0);
    float tot = 0.0;
    float a = 1.0;
    for( int i=0; i<4; i++)
    {
        waterPos += flowOffset * WATER_SPEED;
        flowOffset *= -0.75;
        vec3 v = noised( waterPos ).yzx; // v.xy = deriv, v.z = value
        f += v * a;
        waterPos += v.xy * gradAscent;
        waterPos *= 2.0;
        tot += a;
        a *= foamScale;
    }
    return f / tot;
}
vec3 waterBaseFlow( vec3 p ) {
    // p = waterPos
    return vec3( 1.0, 0.0, riverCurved(p.x) );
}
float waterFlowDepth( vec3 p ) {
    // depth from waterPos to terrainPos under the water
    // p = waterPos
    return WATER_HEIGHT + sceneMap(p).x - p.y;
}
vec3 waterFlowGradient( vec3 p ) {
    // p = waterPos
    vec3 eps = vec3(0.01, 0.0, 0.0);
    float dx = waterFlowDepth( p + eps.xyy ) - waterFlowDepth( p - eps.xyy );
    float dz = waterFlowDepth( p + eps.yyx ) - waterFlowDepth( p - eps.yyx );
    return vec3( dx, 0.0, dz );
}
vec3 waterFlowRate( vec3 p )
{
    // water flow...
    vec3 baseFlow = waterBaseFlow( p );
    vec3 flow = baseFlow;

	float depth = waterFlowDepth( p );
    vec3 dFlow = waterFlowGradient( p );
    
    flow += -dFlow * 40.0 / (1.0 + depth * 1.5);
    flow *= 1.0 / (1.0 + depth * 0.5);

#if 1
    float behindObstacle = 0.5 - dot( normalize(dFlow), -normalize(flow)) * 0.5;
    float slowDist = clamp( depth * 5.0, 0.0, 1.0);
    slowDist = mix(slowDist * 0.9 + 0.1, 1.0, behindObstacle * 0.9);
    slowDist = 0.5 + slowDist * 0.5;
    flow *= slowDist;
#endif

    // water foam...
    float foamScale1 = 0.5;//0.0 ~ 1.0 (default: 0.5)
    float foamScale2 = 0.35;
    float foamCutoff = 0.4;

    float foam = abs(length( flow.xz )) * foamScale1;
	foam += saturate( foam - foamCutoff );
    foam = 1.0 - pow( depth, foam * foamScale2 );
    //foam = 0.1 * foam / depth; // force foam intensity to increase...
    return vec3( flow.xz * 0.6, foam  );
}
vec4 waterNormal( vec3 waterPos, vec3 flowOffset, float foam )
{
    vec2 dWaterPos = max(abs(dFdx(waterPos.xz)), abs(dFdy(waterPos.xz)));
  	float flowScale = max(dWaterPos.x, dWaterPos.y);
    flowScale = (1.0 / (1.0 + flowScale * flowScale * 2000.0));

    float gradAscent = 0.25 - (foam * 1.5);
    vec3 dxy = waterNoise(waterPos.xz * 20.0, flowOffset.xz * 20.0, (0.75 + foam*0.25), gradAscent);
    flowScale *= max(0.25, 1.0 - foam * 5.0); // flatten normal in foam
    vec3 blended = mix( vec3(0.0, 1.0, 0.0), normalize( vec3(dxy.x, flowOffset.y, dxy.y) ), flowScale );
    return vec4( normalize( blended ), dxy.z * flowScale );
}
float waterFoam( vec3 waterPos, vec3 flowOffset, float foam )
{
    // foam = not used...
    float f = waterNoise(waterPos.xz*30.0, flowOffset.xz*50.0, 0.8, -0.5 ).z;
    float amount = 0.2;
    f = max( 0.0, (f - amount) / amount );
    return pow( 0.5, f );
}
vec4 waterFlowingNormal( vec3 waterPos, vec3 flowRate, float foam, float time, out float flowFoam )
{
    float fMag = 2.5 / (1.0 + dot( flowRate, flowRate ) * 5.0);
    float t0 = fract( time );
    float t1 = fract( time + 0.5 );

    float o0 = t0 - 0.5;
    float o1 = t1 - 0.5;
    float weight = abs( t0 - 0.5 ) * 2.0;

    vec3 flowOffset0 = vec3(flowRate.x*o0, fMag, flowRate.z*o0);
    vec3 flowOffset1 = vec3(flowRate.x*o1, fMag, flowRate.z*o1);
    vec4 normal0 = waterNormal( waterPos, flowOffset0, foam );
    vec4 normal1 = waterNormal( waterPos, flowOffset1, foam );
    vec4 normal = mix( normal0, normal1, weight );
    normal.xyz = normalize(normal.xyz);

    o0 *= 0.25;
    o1 *= 0.25;
    flowOffset0 = vec3(flowRate.x*o0, 0.0, flowRate.z*o0);
    flowOffset1 = vec3(flowRate.x*o1, 0.0, flowRate.z*o1);
    float foam0 = waterFoam( waterPos, flowOffset0, foam );
    float foam1 = waterFoam( waterPos, flowOffset1, foam );
    flowFoam = mix( foam0, foam1, weight );

    return normal;
}
vec3 waterEnvColor( vec3 sky, vec3 rd, float gloss )
{
    // WATER_GLOSS_COLOR : appears strongly when view is parallel to the water surface
    const vec3 WATER_GLOSS_COLOR = vec3(0.3, 0.2, 0.2);
    return mix( WATER_GLOSS_COLOR, sky*4.0, saturate(rd.y * (1.0 - gloss*0.5)*0.5 + 0.5) );
}
vec4 waterColor( in vec3 lc, in vec3 ld, in vec3 sky, in vec3 ro, in vec3 rd, in float sceneDist, in float FAR )
// sky = sky color
// sceneDist = rayMarching().x (which used to check if water is visible)
// return xyz = (waterColor), w = (dist from ro to waterPos)
// if water is invisible ==> return vec4(0,0,0, dist to waterPos)
{
    float t = (WATER_HEIGHT-ro.y) / rd.y;
    t = (t > 0.0)? t : FAR;

    vec3 p = ro + rd * t; // waterPos
    gl_FragDepth = getFragDepth( p );

    // flowNormal
    vec3 flowRateFoam = waterFlowRate( p );
    vec3 flowRate = vec3(flowRateFoam.x, 0.0, flowRateFoam.y);
    float flowFoam;
    float foamScale = 1.5;
    float foamOffset = 0.2;
    float foam = saturate( (flowRateFoam.z - foamOffset) * foamScale );
    foam = foam * foam * 0.5;
    vec4 flowNormal = waterFlowingNormal( p, flowRate, foam, time, flowFoam );
    
    if( rd.y < -0.01 )
    {
        t -= (0.04 * (1.0 - flowNormal.w) / rd.y);
    } // here, t = dist from ro to waterPos

    // water visible
    if( t >= sceneDist ) return vec4(0.0, 0.0, 0.0, t);

    // water material
    vec3 albedo = vec3(1.0);
    vec3 specF0 = vec3(0.0204); // F0(water) = 0.0204
    vec2 dp = max(abs(dFdx(p.xz)), abs(dFdy(p.xz)));
    float gloss = max(dp.x, dp.y);
    gloss = exp2( -gloss * 0.3 );

    // flowNormal ==> waterNormal
    vec3 n = normalize( flowNormal.xyz + ld * foam ); // would rather have SSS for foam

    vec3 diffuseCol = vec3(0.0);
    vec3 specularCol = vec3(0.0);

    // waterSpecular ==> specularCol
    vec3 waterDiffuse;
    vec3 waterSpecular;
    getCookShading( albedo, gloss, lc, ld, n, rd, waterDiffuse, waterSpecular );
    specularCol += waterSpecular;

    // reflectCol ===> specularCol
    vec3 reflectCol = reflectRayColor( lc, ld, rd, p, n, FAR ).xyz;
    vec3 reflectRd = reflect( rd, n );
    reflectCol = mix( waterEnvColor(sky, reflectRd, gloss), reflectCol, pow(gloss, 40.0) );
    specularCol += reflectCol;

    // transmitCol + waterDiffuse ==> diffuseCol
    vec3 transmitCol = refractRayColor( lc, ld, rd, p, n, 1.0 / 1.3333, WATER_OPACITY*6.0, FAR ).xyz;
    float foamBlend = 1.0 - pow(flowFoam, foam*5.0);
    diffuseCol = mix(transmitCol, waterDiffuse*0.8, foamBlend );

    // fresnel
    vec3 fresnel = fresnelGloss( n, -rd, specF0, gloss );
    float specScale = saturate(1.0 - foamBlend*4.0);

    // diffuseCol + specularCol ===> result
    vec3 result = mix( diffuseCol, specularCol, fresnel*specScale );
    return vec4(result, t);
}
#endif

//==============================================================================
// seaColor() <=== seaFS.glsl
//==============================================================================

#ifdef USE_SEA_1

const float SEA_CHOPPY = 4.0;
const float SEA_SPEED = 0.8;
const float SEA_FREQ = 0.16;
const float SEA_HEIGHT = 0.6;
const vec3 SEA_BASE_COLOR = vec3(0.1, 0.19, 0.22);
const vec3 SEA_WATER_COLOR = vec3(0.8, 0.9, 0.6);
const mat2 SEA_ROTM2 = mat2(1.6, 1.2, -1.2, 1.6);

float seaOctave( vec2 uv, float choppy )
{
    uv += (2.0*noise(uv)-1.0);
    vec2 wv = 1.0 - abs( sin(uv) );
    vec2 swv = abs( cos(uv) );
    wv = mix( wv, swv, wv );
    return pow( 1.0 - pow(wv.x * wv.y, 0.65), choppy );
}
float seaMap( vec3 p )
{
    float seaTime = time * SEA_SPEED;
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    float d, h = 0.0;
    for(int i = 0; i < 3; i++)
    {
        d = seaOctave(( uv + seaTime)*freq, choppy );
        d += seaOctave( (uv - seaTime)*freq, choppy );
        h += d * amp;
        uv *= SEA_ROTM2;
        freq *= 1.9; amp *= 0.22;
        choppy = mix( choppy, 1.0, 0.2 );
    }
    return p.y - h;
}
float seaMapH( vec3 p )
{
    float seaTime = time * SEA_SPEED;
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    float d, h = 0.0;
    for(int i = 0; i < 5; i++)
    {
        d = seaOctave( (uv + seaTime)*freq, choppy );
        d += seaOctave( (uv - seaTime)*freq, choppy );
        h += d * amp;
        uv *= SEA_ROTM2;
        freq *= 1.9; amp *= 0.22;
        choppy = mix( choppy, 1.0, 0.2 );
    }
    return p.y - h;
}
vec3 seaNormal( vec3 p, float eps )
{
    vec3 n;
    n.y = seaMapH( p );
    n.x = seaMapH( vec3(p.x+eps, p.y, p.z) )     - n.y;
    n.z = seaMapH( vec3(p.x,     p.y, p.z+eps) ) - n.y;
    n.y = eps;
    return normalize(n);
}
float seaTracing( in vec3 ro, in vec3 rd, in float tmax )
{
    float tm = 0.0;
    float tx = tmax;//1000.0;
    float hx = seaMap(ro + rd * tx);
    if( hx > 0.0 ) return tx;
    float hm = seaMap(ro + rd * tm);
    float tmid = 0.0;
    for(int i = 0; i < 8; i++)
    {
        tmid = mix(tm, tx, hm/(hm-hx));
        vec3 p = ro + rd * tmid;
    	float hmid = seaMap(p);
		if( hmid < 0.0 )
        {
        	tx = tmid;
            hx = hmid;
        }
        else
        {
            tm = tmid;
            hm = hmid;
        }
    }
    return tmid;
}
float seaLightDiffuse( vec3 n, vec3 ld, float power ) {
    return pow( dot(n, ld)*0.4 + 0.6, power );
}
float seaLightSpecular( vec3 n, vec3 ld, vec3 rd, float s ) {
    float nrm = (s + 8.0) / (PI * 8.0);
    return pow( max( dot( reflect(rd, n), ld ), 0.0 ), s ) * nrm;
}
vec3 seaColor( vec3 p, vec3 n, vec3 ld, vec3 ro, vec3 rd )
// p = sea position
// n = normal at p
// ld = light direction
// ro, rd = ray definition
{
    float fresnel = 1.0 - max(dot(n,-rd), 0.0);
    fresnel = pow(fresnel, 3.0) * 0.65;

    vec3 reflected = skyColor( ld, reflect(rd, n), 1.0 );
    vec3 refracted = SEA_BASE_COLOR + seaLightDiffuse(n, ld, 80.0) * SEA_WATER_COLOR * 0.12;
    vec3 color = mix( refracted, reflected, fresnel );

    vec3 dist = p - ro;
    float atten = max( 1.0 - dot(dist, dist) * 0.002, 0.0 );//0.001
    color += SEA_WATER_COLOR * (p.y - SEA_HEIGHT) * 0.18 * atten;

    color += vec3( seaLightSpecular( n, ld, rd, 60.0 ) );

    color = pow(color, vec3(1.65)); // we will apply LinearToGamma(2.2) as post-processing...
    return color;
}
vec3 seaColor( vec3 ld, vec3 ro, vec3 rd, float tmax )
{
    float t = seaTracing( ro, rd, tmax );
    vec3 p = ro + rd*t;
    vec3 dist = p - ro;
    float distpp = dot(dist,dist)*0.1 / resolution.x; // estimate "distance per pixel"
    vec3 n = seaNormal( p, distpp );
    //
    gl_FragDepth = getFragDepth( p );
    //
    return seaColor( p, n, ld, ro, rd );
}
#endif

//==============================================================================
// applyDustWind() <=== terrainFS.glsl + cavesFS.glsl
//==============================================================================

#ifdef USE_DUSTWIND_1
float dustWindMap( in vec3 p, in float d, float height, float wind )
{
    p.x += time;
    p.z += time*0.5;
    return triNoise13( p*wind/(d+1.0) ) * smoothstep( height, 0.0, p.y );
}
void applyDustWind( inout vec3 col, in vec3 ro, in vec3 rd, in float t, float amount, float height, float wind )
// t      = distance from ro to hitted position
// amount = dust amount (0.0 ~ 1.0)
// height = dust moves between 0 and height
// wind   = wind turbulency (laminar: 0.01 ~ 0.2, turbulent: 0.3 ~ 1.0)
{
    vec3 dust = vec3(0.85, 0.65, 0.5);
    float d = 0.5;
    for(int i = 0; i < 7; i++)
    {
        vec3 p = ro + rd*d;
        float w = dustWindMap(p, d, height, wind); // w = dust intensity
        col = mix( col, dust, amount*saturate( w * smoothstep(d, d*1.8, t)) );
        d *= 1.8;
        if( d > t ) break;
    }
}
#endif

//==============================================================================
// terrainColor() <=== terrainFS.glsl
//==============================================================================

#ifdef USE_TERRAIN_1

const mat2 rotMat2 = mat2(1.6,-1.2,1.2,1.6);
const float TERRAIN_SCALE = 0.1;  // terrain scene scale: 0.075(high/snow), 0.1, 0.2, 0.3(lower/green)
const float TERRAIN_FREQ = 0.03;  // terrain noise frequency: 0.02, 0.03, 0.04
const float SEA_LEVEL = -2.0;     // -5.0 ~ 0.0

float terrainH( in sampler2D tex, in vec2 x )
{
    vec2  p = x*TERRAIN_FREQ;
    float a = 0.0;
    float b = 1.0;
	vec2  d = vec2(0.0);
    //for(int i = 0; i < 15; i++)
    for(int i = 0; i < 13; i++)
    {
        vec3 n = noised(tex, p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
		b *= 0.5;
        p = rotMat2*p;
    }
    return a/TERRAIN_SCALE + SEA_LEVEL;
}
float terrainM( in sampler2D tex, in vec2 x )
{
    vec2  p = x*TERRAIN_FREQ;
    float a = 0.0;
    float b = 1.0;
	vec2  d = vec2(0.0);
    for(int i = 0; i < 9; i++)
    {
        vec3 n = noised(tex, p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
		b *= 0.5;
        p = rotMat2*p;
    }
    return a/TERRAIN_SCALE + SEA_LEVEL;
}
float terrainL( in sampler2D tex, in vec2 x )
{
    vec2  p = x*TERRAIN_FREQ;
    float a = 0.0;
    float b = 1.0;
	vec2  d = vec2(0.0);
    for(int i = 0; i < 2; i++)//3
    {
        vec3 n = noised(tex, p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
		b *= 0.5;
        p = rotMat2*p;
    }
    return a/TERRAIN_SCALE + SEA_LEVEL;
}
float terrainMap( in sampler2D tex, in vec3 p )
{
    return p.y - terrainM( tex, p.xz );
}
vec3 terrainNormal( in sampler2D tex, in vec3 p, in float t )
{
    // terrainMap:    h(p) = p.y - terrainM( p.xz )
    // terrainNormal: grad h(p) = (dh/dp.x, dh/dp.y, dh/dp.z)
    vec2 eps = vec2( 0.002*t, 0.0 );
    return normalize( vec3( terrainH(tex, p.xz-eps.xy) - terrainH(tex, p.xz+eps.xy),
                            2.0*eps.x,
                            terrainH(tex, p.xz-eps.yx) - terrainH(tex, p.xz+eps.yx) ) );
}
vec3 terrainColor( in sampler2D tex, in vec3 lc, in vec3 ld, in vec3 p, in float t, in float tmax )
// tex = 'images/raymarch/grayNoise256.png'
// p = terrainPos
// t = dist from ro(cameraPos) to p
// tmax = FAR
{
    vec3 n = terrainNormal( tex, p, t );

    float r = texture( tex, (7.0/TERRAIN_SCALE)*p.xz/256.0 ).x;

    // soil
    vec3 soilA = vec3(0.08,0.05,0.03);
    vec3 soilB = vec3(0.10,0.09,0.08);
    vec3 col = (0.75 + 0.25*r)*mix( soilA, soilB, texture(tex,0.0007*vec2(p.x,p.y*19.19)/TERRAIN_SCALE).x );
    // rock
    vec3 rock = 0.2*vec3(0.45, 0.30, 0.15);
    col = mix( col, rock*(0.5 + 0.5*r), smoothstep(0.85, 0.9, n.y) );
    // weed
    vec3 weed = 0.1*vec3(0.30, 0.30, 0.10);
    col = mix( col, weed*(0.5 + 0.5*r), smoothstep(0.9, 0.95, n.y) );
    // grass
    vec3 grass = 0.35*vec3(0.16, 0.3, 0.0);
    col = mix( col, grass*(0.25 + 0.75*r), smoothstep(0.95, 1.0, n.y) );
    // snow
    float hmin = 10.0/TERRAIN_SCALE;
    float hmax = 30.0/TERRAIN_SCALE;
    applySnow( col, tex, hmin, hmax, p/TERRAIN_SCALE, n );

    // lighting
    float amb = saturate( 0.5 + 0.5*n.y );
    float dif = saturate( dot( ld, n ) );
    float bac = saturate( 0.2 + 0.8*dot(normalize(vec3(-ld.x, 0.0, ld.z)), n) );
    float shd = (dif >= 0.001)? sceneShadow(p + ld*0.01/TERRAIN_SCALE, ld, 0.01/TERRAIN_SCALE, tmax/TERRAIN_SCALE, 2.0) : 1.0;
    vec3 cshd = pow( vec3(shd), vec3(1.0, 1.2, 1.5) ); // colorize shadow penumbras
    vec3 lin = dif*vec3(7.0,5.0,3.0)*1.3*cshd * lc;
    lin += amb*vec3(0.4,0.6,1.0)*1.2;
    lin += bac*vec3(0.4,0.5,0.6);
    col *= lin;

    return col;
}
#endif

//==============================================================================
// terrainColor() <=== forestFS.glsl
//==============================================================================

#ifdef USE_TERRAIN_2

float terrainNoise( in sampler2D tex, in vec2 x )
{
    float f = 1.92;//2.0;
    float a = 0.0;
    float b = 0.5;
    for(int i=0; i<6; i++)
    {
        float n = noise(tex, x);
        a += b * n;
        b *= 0.55;//0.55;
        x = f * rotM2 * x;
    }
	return a;
}
float terrainH( in sampler2D tex, in vec2 p )
{
    const float terrain_freq = 0.01;
    const float terrain_scale = 70.0;//10.0(flat) ~ 100.0(high)
    p *= terrain_freq;
    float e = -1.0 + 2.0*terrainNoise( tex, p + vec2(1.0,-2.0) );
    e *= terrain_scale;
    return e;
}
float terrainM( in sampler2D tex, in vec2 p )
{
    return terrainH( tex, p );
}
float terrainL( in sampler2D tex, in vec2 p )
{
    return terrainH( tex, p );
}

float terrainMap( in sampler2D tex, in vec3 p )
{
    return p.y - terrainM( tex, p.xz );
}
vec3 terrainNormal( in sampler2D tex, in vec3 p, in float t )
{
    vec2 eps = vec2(0.002*t, 0.0);
	return normalize( vec3(terrainH(tex, p.xz-eps.xy) - terrainH(tex, p.xz+eps.xy),
                           2.0*eps.x,
                           terrainH(tex, p.xz-eps.yx) - terrainH(tex, p.xz+eps.yx) ) );
}
float terrainShadow( in sampler2D tex, in vec3 ro, in vec3 rd, in float tmin )
{
    float shade = 1.0;
    float t = tmin;
    for(int i=0; i<8; i++)//32
    {
        vec3  pos = ro + t*rd;
        float h = terrainMap( tex, pos );
        shade = min( shade, 32.0*h/t );
        if( shade < 0.0001 ) break;
        t += clamp( h, 1.0+t*0.1, 50.0 );
    }
    return saturate( shade );
}
vec3 terrainColor( in sampler2D tex, in vec3 ld, in vec3 rd, in vec3 p, in float t )
// tex = noise texture (eg: grayNoise256.png)
// ld = light direction (eg: SUN_LIGHT)
{
    vec3 n = terrainNormal( tex, p, t );
    // bump map
    #if 1
        n = normalize( n + 1.28*(1.0-abs(n.y)) * fbmd_7(p*0.3*vec3(1.0,0.2,1.0)).yzw );//1.28=0.8*0.8*2.0
    #endif

    vec3 col = vec3(0.18,0.11,0.10)*0.75;
    col = mix( col, vec3(0.1,0.1,0.0)*0.3, smoothstep(0.7, 0.9, n.y) );
    #if 0
        col *= 1.0 + 2.0*fbm_4( iChannel0, p*0.2*vec3(1.0,4.0,1.0) );
    #endif
    
    float sha = 0.0;
    float dif = saturate( dot(n, ld) );
    if( dif > 0.0001 ) 
    {
        sha = terrainShadow( tex, p+n*0.01, ld, 0.01 );
        dif *= sha;
    }
    vec3  ref = reflect(rd, n);
    float bac = saturate( dot(normalize(vec3(-ld.x, 0.0, -ld.z)), n) ) * saturate( (p.y+100.0)/100.0 );
    float dom = saturate( 0.5 + 0.5*n.y );
    vec3  lin  = 0.2*mix(0.1*vec3(0.1,0.2,0.0), vec3(0.7,0.9,1.0), dom);//pow(vec3(occ),vec3(1.5,0.7,0.5));
    lin += 5.0*vec3(1.0,0.9,0.8)*dif;        
    lin += 0.35*vec3(1.0)*bac;
    col *= lin;

    #if 1
        applyFog( col, FOG_COLOR, 0.0005, t );
    #endif

    return col;
}
#endif

//==============================================================================
// terrainColor() <=== canyonFS.glsl
//==============================================================================

#ifdef USE_TERRAIN_3

float terrainDetails( in sampler2D detailTex, in vec3 p )
{
#if 0
    float f;
    f  = 0.5000*noise( detailTex, p ); p = rotM3*p*2.02;
    f += 0.2500*noise( detailTex, p ); p = rotM3*p*2.03;
    f += 0.1250*noise( detailTex, p ); p = rotM3*p*2.01;
    f += 0.0625*noise( detailTex, p );
    return f;
#else
	return fbm_4( detailTex, p );
#endif
}
float terrainH( in sampler2D heightTex0, in sampler2D heightTex1, in vec2 q )
// heightTex0: (global) height texture
// heightTex1:  (local) height texture
// shapeNoise: 0.1=(almost_flat), 0.5=(smooth+wide), 2.5=(sharp+rugged)
// seaLevel: height from the bottom of terrain
// riseHeight (-15.0 ~ 15.0): rise above seaLevel (if negative, terrain exists below seaLevel)
// sinkHeight (0.0 ~ 1.0): depth below seaLevel
{
	// shapeNoise (0.1 ~ 2.5)
	float shapeNoise = 1.0;//1.5;//1.0;
	q *= shapeNoise;

	float th = smoothstep( 0.0, 0.7, textureLod( heightTex0, 0.001*q, 0.0 ).x );    // heightTex0 = 1024 x 1024
    float rr = smoothstep( 0.1, 0.5, textureLod( heightTex1, 2.0*0.03*q, 0.0 ).y ); // heightTex1 = 1024 x 1024

	// seaLevel
	float seaLevel = 0.5;
	float h = 0.7 - seaLevel;

	h -= -0.15 + (1.0-0.6*rr)*(1.5-1.0*th) * 0.3*(1.0-textureLod( heightTex0, 0.04*q*vec2(1.2,0.5), 0.0 ).x);
	//h += -0.15 + (1.0-0.6*rr)*(1.5-1.0*th) * 0.3*(1.0-textureLod( heightTex0, 0.04*q*vec2(1.2,0.5), 0.0 ).x);

	// riseHeight (-15.0 ~ 15.0)
	float riseHeight = 7.0;
	h += th * riseHeight;

	// sinkHeight (0.0 ~ 1.0)
	float sinkHeight = 0.3;
	h -= rr * sinkHeight;
    return h;
}
float terrainL( in sampler2D heightTex0, in vec2 q )
{
	float th = smoothstep( 0.0, 0.7, textureLod( heightTex0, 0.001*q, 0.0 ).x );
	float h = 0.2;//1.0 0.2
	h += th * 7.0;//7.0
	return h;
}
float terrainMap( in sampler2D heightTex0, in sampler2D heightTex1, in sampler2D detailTex, in vec3 p )
{
	// base height
	float h = terrainH( heightTex0, heightTex1, p.xz );

	// added: details
	float detailAmount = 0.001;//0.15;//0.25;	// 0.0 ~ 0.5
    float layerAmount = 1.5;//3.0;		// 0.0 ~ 3.0
    float displaceAmount = 1.5;//3.0;	// 0.0 ~ 5.0

    float d = terrainDetails( detailTex, detailAmount * p * vec3(1.0, layerAmount, 1.0) );
    d *= displaceAmount;

	return (p.y - h + d) * 0.25;
}
vec3 terrainNormal( in sampler2D heightTex0, in sampler2D heightTex1, in vec3 p, in float t )
{
	vec2 eps = vec2( 0.002*t, 0.0 );
    return normalize( vec3( terrainH(heightTex0, heightTex1, p.xz-eps.xy) - terrainH(heightTex0, heightTex1, p.xz+eps.xy),
                            2.0*eps.x,
                            terrainH(heightTex0, heightTex1, p.xz-eps.yx) - terrainH(heightTex0, heightTex1, p.xz+eps.yx) ) );
}
float terrainShadow( in sampler2D heightTex0, in sampler2D heightTex1, in sampler2D detailTex, in vec3 ro, in vec3 rd, float tmin )
{
	float shade = 1.0;
    float t = tmin;
    for(int i=0; i<32; i++)//8 64
    {
        float h = terrainMap(heightTex0, heightTex1, detailTex, ro + rd * t);
        shade = min( shade, 32.0*h/t );//32.0
        t += clamp( h, 0.5, 1.0 );
		if( h < 0.001 ) break;
    }
    return min( max(shade, 0.0) + 0.05, 1.0 ); // 0.3 or 0.1 (preference)
}
float terrainAO( in sampler2D heightTex0, in sampler2D heightTex1, in sampler2D detailTex, in vec3 p, in vec3 n )
{
	float ao = 0.0;
    float s = 1.0;
    for(int i=0; i<2; i++)//6
    {
        float off = 0.001 + 0.2 * float(i)/5.0;
        float t = terrainMap( heightTex0, heightTex1, detailTex, n * off + p );
        ao += ( off - t ) * s;
        s *= 0.4;
    }
    return smoothstep( 0.0, 1.0, clamp(1.0-12.0*ao, 0.0, 1.0) );
}
vec3 terrainColor( in sampler2D heightTex0, in sampler2D heightTex1, in sampler2D detailTex, in vec3 ld, in vec3 rd, in vec3 p, in float t )
{
	vec3 n = terrainNormal( heightTex0, heightTex1, p, t );

	float gloss = 2.0;//1.0; // 0.1(trash-heap) ~ 10.0(smooth)
	vec3 uvw = gloss * p;

	// bump normal
	n = getBumpNormal( heightTex0, uvw, n, 0.075 );

	// materials
	vec3 te = 0.05 + texCube( heightTex0, 0.15*uvw, n ).xyz;
	vec4 mate;
	mate.xyz = 0.6*te;				// material diffuse
	mate.w = 1.5*(0.5+0.5*te.x);	// material shininess

	// added: soil color using heightTex0
	float th = smoothstep( 0.1, 0.4, texCube( heightTex0, 0.002*uvw, n ).x );
	vec3 dcol = mix( vec3(0.2, 0.3, 0.0), 0.2*vec3(0.65, 0.4, 0.2), 0.2+0.8*th );
	mate.xyz = mix( mate.xyz, 2.0*dcol, th*smoothstep(0.0, 1.0, n.y) );

	// added: snow(white) using heightTex1
	float rr = smoothstep( 0.2, 0.4, texCube( heightTex1, 0.04*uvw, n ).y );
	mate.xyz *= mix( vec3(1.0), 2.25*vec3(0.25,0.24,0.22), rr );
	mate.xyz *= 1.5*pow(texCube( heightTex1, 8.0*uvw, n ).xyz, vec3(0.5));
	mate = mix( mate, vec4(vec3(1.0), 0.0), smoothstep(0.8, 0.9, n.y + n.x*0.6*te.x*te.x) );
	mate.xyz *= 1.5;

	// lighting
	float sky = 0.0;
	sky += 0.2*orenNayarDiffuse( normalize(vec3( 0.0, 1.0, 0.0 )), n, -rd, 1.0 );
	sky += 0.2*orenNayarDiffuse( normalize(vec3( 3.0, 1.0, 0.0 )), n, -rd, 1.0 );
	sky += 0.2*orenNayarDiffuse( normalize(vec3(-3.0, 1.0, 0.0 )), n, -rd, 1.0 );
	sky += 0.2*orenNayarDiffuse( normalize(vec3( 0.0, 1.0, 3.0 )), n, -rd, 1.0 );
	sky += 0.2*orenNayarDiffuse( normalize(vec3( 0.0, 1.0,-3.0 )), n, -rd, 1.0 );
	float dif = orenNayarDiffuse( ld, n, -rd, 1.0 );
	vec3 blig = normalize(vec3(-ld.x, 0.0, -ld.z));
	float bac = orenNayarDiffuse( blig, n, -rd, 1.0 );

	float sha = 0.0;
	if( dif > 0.001 ) sha = terrainShadow( heightTex0, heightTex1, detailTex, p+0.01*n, ld, 0.005 );
	float spe = mate.w * pow(saturate(dot(reflect(rd,n),ld)), 2.0) * saturate(dot(n,ld));
	float occ = terrainAO( heightTex0, heightTex1, detailTex, p, n );
	vec3 lin = vec3(0.0);
	lin += 7.0*dif*vec3(1.20,0.50,0.25)*vec3(sha,sha*0.5+0.5*sha*sha, sha*sha);
	lin += 1.0*sky*vec3(0.10,0.50,0.70)*occ;
	lin += 2.0*bac*vec3(0.30,0.15,0.15)*occ;
	lin += 0.5*vec3(spe)*sha*occ;
	vec3 col = mate.xyz * lin;
	return col;
}
#endif

//==============================================================================
// treesColor() <=== forestFS.glsl
//==============================================================================

#ifdef USE_TREES_1

const float TREES_HEIGHT = 2.0;

float terrainM( in sampler2D tex, in vec2 p );
vec3 terrainNormal( in sampler2D tex, in vec3 p, in float t );

float treesMap( in sampler2D tex, in vec3 p )
// tex = noise texture (eg: grayNoise256.png)
{
    //float base = terrainM(p.xz);
    float base = terrainM(tex, p.xz)*0.925;

    float d = 20.0;//10.0
    vec2 n = floor( p.xz );
    vec2 f = fract( p.xz );
    for(int j=-1; j<=1; j++)
    for(int i=-1; i<=1; i++)
    {
        vec2  g = vec2( float(i), float(j) );
        vec2  o = hash22( n + g );
        vec2  v = hash22( n + g + vec2(13.1,71.7) );
        vec2  r = g - f + o;

        float height = TREES_HEIGHT * (0.4 + 0.8*v.x);
        float width = 0.9*(0.5 + 0.2*v.x + 0.3*v.y);
        vec3  q = vec3(r.x, p.y-base-height*0.5, r.y);
        #if 1
            float k = fEllipsoid( q, vec2(width,0.5*height).xyx );
        #else
            vec3 a = vec3(0.0, -height*0.5, 0.0);
            vec3 b = vec3(0.0, height*0.5, 0.0);
            float k = fCapsule( q, a, b, width );
        #endif
        d = min( d, k );
    }

    // distort ellipsoids to make them look like trees (works only in the distance really)
    float rt = 350.0 * rand( p.xz );
    //if( rt < 350.0 )
    {
        float s = -1.0 + 2.0*fbm_4( tex, p*3.0 );
        float att = 1.0-smoothstep(150.0, 350.0, rt);
        d += 2.0*s*s*att*att;
    }
    
    return d;
}
vec3 treesNormal( in sampler2D tex, in vec3 p, in float t )
{
    vec2 e = vec2(t,-t) * 0.002;
    return normalize( e.xyy * treesMap(tex, p + e.xyy) 
                    + e.yyx * treesMap(tex, p + e.yyx) 
                    + e.yxy * treesMap(tex, p + e.yxy) 
                    + e.xxx * treesMap(tex, p + e.xxx) );
}
float treesShadow( in sampler2D tex, in vec3 ro, in vec3 rd )
{
    float shade = 1.0;
    float t = 0.02;
    //for( int i=0; i<50; i++ )
    for( int i=0; i<10; i++ )
    {
        float h = treesMap( tex, ro + rd*t );
        shade = min( shade, 32.0*h/t );
        t += h;
        if( shade < 0.001 || t > 20.0 ) break;
    }
    return saturate( shade );
}
vec3 treesColor( in sampler2D tex, in vec3 ld, in vec3 pos, in vec3 rd, float t )
// tex = noise texture (eg: grayNoise256.png)
// ld = light direction (eg: SUN_LIGHT)
{
    vec3 terrainN = terrainNormal( tex, pos, t );

    vec3 treesN = treesNormal( tex, pos, t );
    vec3 n = normalize( treesN + 2.5*terrainN );

    // lighting
    float sha = 1.0;
    vec3  ref = reflect(rd, n);
    float occ = 1.0;
    float dif = saturate(0.1 + 0.9*dot(n, ld));
    if( dif > 0.0001 )
    {
        sha *= treesShadow( tex, pos+n*0.1, ld ); // only cast in non-terrain-occluded areas
    }
    float dom = saturate( 0.5 + 0.5*n.y );
    float fre = saturate( 1.0 + dot(n,rd) );
    float spe = pow(saturate(dot(ref, ld)), 9.0) * dif*sha * (0.2 + 0.8*pow(fre, 5.0)) * occ;

    // lights
    vec3 lin = 0.5*mix(0.1*vec3(0.1,0.2,0.0),vec3(0.6,1.0,1.0),dom*occ);
    lin += 10.0*vec3(1.0,0.9,0.8)*dif*occ*sha;
    lin += 0.5*vec3(0.9,1.0,0.8)*pow(fre, 3.0)*occ;
    lin += 0.05*vec3(0.15,0.4,0.1)*occ;
   
    // material
    float brownAreas = fbm_4( tex, pos.zx*0.03 );
    vec3 col = vec3(0.08,0.09,0.02);
    col = mix( col, vec3(0.06,0.05,0.01)*1.1, 1.0-smoothstep(0.9,0.91,terrainN.y) );
    col = mix( col, vec3(0.25,0.16,0.01)*0.15, 0.7*smoothstep(0.1,0.3,brownAreas)*smoothstep(0.5,0.8,terrainN.y) );
    col *= 1.6;

    // brdf * material
    col *= lin;
    col += spe*1.2*vec3(1.0,1.1,2.5);

    #if 1
        applyFog( col, FOG_COLOR, 0.0005, t );
    #endif

    return col;
}
#endif

#endif // RAYMARCH_SCENE