@import ./include;
@import ./utils;
@import ./noises;
@import ./distances;
#define USE_DUSTWIND_1
@import ./scene;

#define map( pos )  sceneMap( pos ).x

// textureMaps[0] = 'images/raymarch/pebbles.jpg'
// textureMaps[1] = 'images/raymarch/organic1.jpg'
#define iChannel0   textureMaps[0]
#define iChannel1   textureMaps[1]

//==============================================================================
// Desert Canyon
// The far plane. I'd like this to be larger, but the extra iterations required to render the 
// additional scenery starts to slow things down on my slower machine.

// Choose...
#define ENABLE_AUTO_VIEW
#define ENABLE_DUST_WIND
#define ENABLE_PHONG_LIGHTING

// Choose...
//#define DESERT_CAVE
//#define SMOOTH_CAVE
//#define HELL_CAVE
#define CAVE_CANYON         // similar to DESERT_CAVE (default)
//#define FLOATING_CANYON

#define FAR 65.0

// used in path function & used to shape the tunnel and guide the camera.
const float freqA = 0.15/3.75;
const float freqB = 0.25/2.75;
const float ampA = 20.0;
const float ampB = 4.0;

mat2 rot2( float th ){ vec2 a = sin(vec2(1.5707963, 0) + th); return mat2(a, -a.y, a.x); }

// The path is a 2D sinusoid that varies over time, depending upon the frequencies and amplitudes.
vec2 path( in float z )
{
    return vec2( ampA*sin(z * freqA), ampB*cos(z * freqB) + 3.0*(sin(z*0.025) - 1.0) );
}

#ifdef DESERT_CAVE
// The canyon, complete with hills, gorges and tunnels. I would have liked to provide a far
// more interesting scene, but had to keep things simple in order to accommodate slower machines.
float desertCanyon(in vec3 p)
{
    // Indexing into the pebbled texture to provide some rocky surface detatiling. I like this
    // texture but I'd much rather produce my own. From what I hear, Shadertoy will be providing
    // fixed offscreen buffer sizes (like 512 by 512, for instance) at a later date. When that
    // happens, I'll really be able to do some damage. :)

    // rocky surface detailing
    float tx = textureLod(iChannel0, p.xz/16.0 + p.xy/80.0, 0.0).x;
  
    // sinusoidal layers (to produce the rocky hills)
    vec3 q = p*0.25;
    float h = dot(sin(q)*cos(q.yzx), vec3(0.222)) + dot(sin(q*1.5)*cos(q.yzx*1.5), vec3(0.111));

    // The terrain, so to speak. Just a flat XZ plane, at zero height, with some hills added.
    float d = p.y + h*6.0;
  
    // Reusing "h" to provide an undulating base layer on the tunnel walls.
    q = sin(p*0.5 + h);
    h = q.x*q.y*q.z;
  
	// Producing a single winding tunnel. If you're not familiar with the process, this is it.
    // We're also adding some detailing to the walls via "h" and the rocky "tx" value.
    p.xy -= path(p.z);
    float tunnel = 1.5 - length(p.xy*vec2(0.33, 0.66)) + h + (1.0 - tx)*0.25;

	// Smoothly combine the terrain with the tunnel - using a smooth maximum - then add some
    // detailing. I've also added a portion of the tunnel term onto the end, just because
    // I liked the way it looked more. 
    return fsmax(d, tunnel, 2.0) - tx*0.5 + tunnel*0.8;
}
#endif

#ifdef SMOOTH_CAVE
float smoothCave(in vec3 p)
{
    // rocky surface detailing
    float tx = textureLod(iChannel0, p.xz/16.0 + p.xy/80.0, 0.0).x;

    // sinusoidal layers (to produce the rocky hills)
    vec3 q = p*0.25;
    float freq = 1.5;// The larger value, the more empty space in the internal space
    float h = dot(sin(q)*cos(q.yzx), vec3(0.122)) + dot(sin(q*freq)*cos(q.yzx*freq), vec3(0.111));

    // terrain (with some hills added)
    float d = p.y*0.2 + h*20.0;
  
    // undulating base layer on the tunnel walls
    q = sin(p*0.5 + h);
    h = q.x*q.y*q.z;
  
	// winding tunnel (adding details to the walls via "h" and the rocky "tx" value)
    p.xy -= path(p.z);
    float tunnel = 1.5 - length(p.xy*vec2(0.33, 0.66)) + h + (1.0 - tx)*0.25;

	// combine terrain with the tunnel
    //return fsmax(d, tunnel, 2.0) + tx*0.5 + tunnel*0.8;
    return fsmax(d, tunnel, 2.0) + tx*0.5 + tunnel*0.8 + tunnel * tx;
}
#endif

#ifdef HELL_CAVE
float hellCanyon(in vec3 p)
{
    // rocky surface detailing
    float tx = textureLod(iChannel0, p.xz/16.0 + p.xy/80.0, 0.0).x;
  
    // sinusoidal layers (to produce the rocky hills)
    vec3 q = p*0.15;
    float freq = 1.5;// The larger value, the more empty space in the internal space
    float h = dot(sin(q)*cos(q.yzx), vec3(0.122)) + dot(sin(q*freq)*cos(q.yzx*freq), vec3(0.111));

    // terrain (with some hills added)
    float d = p.y*3.1 + h*16.0;

    // undulating base layer on the tunnel walls
    q = sin(p*0.5 + h);
    h = q.x*q.y*q.z;

    // winding tunnel
    p.xy -= path(p.z);
    float tunnel = 1.5 - length(p.xy*vec2(0.33, 0.66)) + h + (1.0 - tx)*0.25;

    // combine terrain with the tunnel
    return fsmax(d, tunnel, 2.0) + tx*0.5 - tunnel*0.8;
}
#endif

#ifdef CAVE_CANYON
float caveCanyon(in vec3 p)
{
    // rocky surface detailing
    float tx = textureLod(iChannel0, p.xz/16.0 + p.xy/80.0, 0.0).x;
  
    // sinusoidal layers (to produce the rocky hills)
    vec3 q = p*0.24;
    float freq = 2.8;
    float h = dot(sin(q)*cos(q.yzx), vec3(0.122)) + dot(sin(q*freq)*cos(q.yzx*freq), vec3(0.111));

    // terrain (with some hills added)
    float d = p.y - h*20.0;

    // undulating base layer on the tunnel walls
    q = sin(p*0.5 + h);
    h = q.x*q.y*q.z;

    // winding tunnel
    p.xy -= path(p.z);
    float tunnel = 1.5 - length(p.xy*vec2(0.33, 0.66)) + h + (1.0 - tx)*0.25;

    // combine terrain with the tunnel
    return fsmax(d, tunnel, 2.0) - tx*0.5 + tunnel*0.8;
}
#endif

#ifdef FLOATING_CANYON
float floatingCanyon(in vec3 p)
{
    // rocky surface detailing
    float tx = textureLod(iChannel0, p.xz/16.0 + p.xy/80.0, 0.0).x;
  
    // sinusoidal layers (to produce the rocky hills)
    vec3 q = p*0.25;
    float freq = 1.5;// The larger value, the more empty space in the internal space
    float h = dot(sin(q)*cos(q.yzx), vec3(0.122)) + dot(sin(q*freq)*cos(q.yzx*freq), vec3(0.111));

    // terrain (with some hills added)
    float d = p.y*0.2 + h*22.0;
  
    // undulating base layer on the tunnel walls
    q = sin(p*0.5 + h);
    h = q.x*q.y*q.z;

    // winding tunnel
    p.xy -= path(p.z);
    float tunnel = 1.5 - length(p.xy*vec2(0.33, 0.66)) + h + (1.0 - tx)*0.25;

    // combine terrain with the tunnel
    return fsmax(d, tunnel, 2.0) + tx*0.5 - tunnel*0.8;
}
#endif

vec2 sceneMap( in vec3 p )
{
#if defined( DESERT_CAVE )
    return vec2( desertCanyon(p), 0.0 );

#elif defined( SMOOTH_CAVE )
    return vec2( smoothCave(p), 0.0 );

#elif defined( HELL_CAVE )
    return vec2( hellCanyon(p), 0.0 );

#elif defined( CAVE_CANYON )
    return vec2( caveCanyon(p), 0.0 );

#elif defined( FLOATING_CANYON )
    return vec2( floatingCanyon(p), 0.0 );

#endif
}

void applyLighting( inout vec3 col, in vec3 p, in vec3 n, in vec3 ld, in vec3 rd )
// ld = light direction from surface to light
{
    float diffuseIntensity = 3.0;
    float specularIntensity = 1.0;//2.25;
    float shininess = 16.0;
    float F0 = 0.2;// F0(dielectric)=0.04, F0(water)=0.0204

    //float shadow = sceneShadow(p, ld, 0.05, FAR, 8.0);
    float shadow = sceneShadow(p+ld*0.05, ld, 0.05, FAR, 8.0);
    float ao = sceneAO(p, n);
    vec3 v = -rd;
    vec3 h = normalize(ld + v);
    float diffuse = max( dot(ld, n), 0.0) * diffuseIntensity;
    float specular = pow( max(dot(n, h), 0.0), shininess) * specularIntensity;
    float F = mix(F0, 1.0, pow(1.0 - max(dot(h,v),0.0), 5.0));// schlick approximation
    float ambient = saturate( 0.3 + 0.1*n.y );
    col = mix( col*(diffuse + ambient)*shadow*ao, vec3(specular*shadow), F );
}

void main()
{
#ifdef ENABLE_AUTO_VIEW
	// Screen coordinates
	vec2 u = (gl_FragCoord.xy - resolution.xy*0.5)/resolution.y;
	
	// Camera Setup.
	vec3 lookAt = vec3(0.0, 0.0, time*8.0);  // "Look At" position.
	vec3 ro = lookAt + vec3(0.0, 0.0, -0.1); // Camera position, doubling as the ray origin.
 
	// Using the Z-value to perturb the XY-plane.
	// Sending the camera and "look at" vectors down the tunnel. The "path" function is 
	// synchronized with the distance function.
	lookAt.xy += path(lookAt.z);
	ro.xy += path(ro.z);

    // Using the above to produce the unit ray-direction vector.
    float FOV = 3.14159/3.; // FOV - Field of view.
    vec3 forward = normalize(lookAt-ro);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);

    // rd - Ray direction.
    vec3 rd = normalize(forward + FOV*u.x*right + FOV*u.y*up);
    
    // Swiveling the camera about the XY-plane (from left to right) when turning corners.
    // Naturally, it's synchronized with the path in some kind of way.
	rd.xy = rot2( path(lookAt.z).x/64. )*rd.xy;
#else
    @import ./ray;
    vec2 u = ndc; // ndc = normalized device coordinates [-1, +1]
#endif

    // light position
    vec3 lp = vec3(FAR*0.5, FAR, FAR) + vec3(0, 0, ro.z);

    // sky & clouds
#ifdef HELL_CAVE
    vec3 sky = vec3(0.0, -0.0, -0.0);
    applyClouds( sky, ro, rd );
#else
    vec3 sky = skyColor(normalize(lp - ro), rd, 1.0);
    applyClouds( sky, ro, rd );
#endif

    vec3 col = sky;

    // ray marching...
    float t = rayTracing( ro, rd, 0.001, FAR, 0.3 );//0.4

    gl_FragDepth = 0.99;
    
    // If we've hit the ground, color it up.
    if( t < FAR )
    {
#ifdef ENABLE_PHONG_LIGHTING
        vec3 p = ro+t*rd;
        vec3 n = sceneNormal( p, t );

        gl_FragDepth = getFragDepth( p );

        // light direction
        vec3 ld = lp - p;
        ld /= max( length(ld), 0.0001 );

        // 0) bump normal (for details)
        const float texScale = 0.2;
        n = getBumpNormal(iChannel1, p*texScale, n, 0.01/(1.0 + t/FAR));

        // 1) base color (soil) (see 'skin peeler' for details)
        col = mix( vec3(0.8, 0.5, 0.3), vec3(0.5, 0.25, 0.125), (p.y + 1.0)*0.15 );
        col = clamp(col, vec3(0.5, 0.25, 0.125), vec3(1.0));

        // 2) texture
        col = smoothstep(-0.5, 1.0, texCube(iChannel1, p*texScale, n).rgb)*(col + 0.25);

        // 3) anisotropic filtering (for crisp textures)
        float crisp = 0.45;
        col = saturate(col + noise13(p*48.0)*crisp - 0.15);

        // 4) color in the crack ==> more darker
        float curv = cheapCurvature(p)*0.9 + 0.1;
        col *= smoothstep(0.0, 0.8, curv);

        // 5) sky reflection ===> almost ineffective
        //col += skyColor(ld, reflect(rd, n), 1.0) * F * 0.2;

        // 6) gamma correction
        col = pow(col, vec3(2.2));

        // 7) lighting
        applyLighting( col, p, n, ld, rd );
#else
        vec3 p = ro+t*rd;
        vec3 n = sceneNormal( p, t );
        gl_FragDepth = getFragDepth( p );

        vec3 ld = lp - p;
        ld /= max( length(ld), 0.0001 );

        const float texScale = 1.0/6.0;
        n = getBumpNormal(iChannel1, p*texScale, n, 0.007/(1.0 + t/FAR));

        float shd = sceneShadow(p, ld, 0.05, FAR, 8.0);
        float ao = sceneAO(p, n);
        float curv = cheapCurvature(p)*0.9 + 0.1;

        vec3 v = -rd;
        vec3 h = normalize(ld + v);
        float dif = max( dot( ld, n ), 0.0);
        float spe = pow( max(dot( n, h ),0.0), 5.0); // 5.0 = shininess
        float fre = saturate(1.0 - dot(v, n)); // (1 - dotNV) for fresnel term
		float schlick = pow(1.0 - max(dot(rd, normalize(rd + ld)), 0.0), 5.0);
		float F = mix(0.2, 1.0, schlick);  // F0 (hard clay) = 0.2
        float amb = fre*F + 0.06*ao;

        col = mix( vec3(0.8, 0.5, 0.3), vec3(0.5, 0.25, 0.125), (p.y + 1.0)*0.15 );
        col = clamp( col, vec3(0.5, 0.25, 0.125), vec3(1.0));
        col = smoothstep(-0.5, 1.0, texCube(iChannel1, p*texScale, n).rgb)*(col + 0.25);
        col = saturate(col + noise13(p*48.0)*0.3 - 0.15);
        col = pow(col, vec3(1.5));
        col *= smoothstep(0.0, 0.7, curv);
        col += skyColor(ld, reflect(rd, n), 1.0)*fre*F*0.5;
        col = (col*(dif + 0.1) + F*spe)*shd*ao + amb*col;
#endif
    }

    // terrain + sky(with fog) ===> foggy terrain
    col = mix(col, sky, sqrt(smoothstep(FAR - 15.0, FAR, t)));

#ifdef ENABLE_DUST_WIND
    #ifdef FLOATING_CANYON
        float dustAmount = 0.005;
    #else
        float dustAmount = 0.01;
    #endif
    float dustHeight = 3.0;
    float windTurbulency = 0.3;
    applyDustWind( col, ro, rd, t, dustAmount, dustHeight, windTurbulency );
#endif

#ifdef ENABLE_PHONG_LIGHTING
    col = pow(col, vec3(0.45)); // LinearToGamma
#else
    col = pow(max(col, 0.0), vec3(0.75));
#endif

    // vignetting
#if 0
    u = gl_FragCoord.xy / resolution.xy;
    col *= pow( 16.0*u.x*u.y*(1.0-u.x)*(1.0-u.y), 0.0625);
#else
    col = Vignetting( col, 0.5 );
#endif

	gl_FragColor = vec4( saturate(col), 1.0 );
}