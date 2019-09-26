@import ./include;

#define texChannel0 textureMaps[0]

// Choose...
#define ENABLE_AUTO_CAMERA

// Choose...
//#define ENABLE_SURF_JAG1      // sharp
#define ENABLE_SURF_JAG2        // default
//#define ENABLE_SURF_JAG3      // flat
//#define ENABLE_SURF_JAG4
//#define ENABLE_SURF_JAG5
//#define ENABLE_SURF_JAG6      // round

// Choose...
#define ENABLE_ROUND_TUNNEL     // default
//#define ENABLE_SQUARE_TUNNEL
//#define ENABLE_ROUNDED_SQUARE_TUNNEL
//#define ENABLE_MINIMALISTS

// Choose...
#define ENABLE_GLOW           // default
//#define ENABLE_SHINY

// Choose...
//#define ENABLE_EDGE_EFFECT

// Choose the camera speed
#define CAMERA_SPEED    5.0
// 1.0 = slow
// 5.0 = normal (default)
// 10.0 = fast

// Choose the ground height
#define GROUND_HEIGHT   1.0
// 0.0 = wide flat groud
// 1.0 = narrow flat ground (default)
// 2.0 = bumpy ground (no flat)
// 3.0 = more bumpy ground


// Grey scale.
float getGrey(vec3 p){ return p.x*0.299 + p.y*0.587 + p.z*0.114; }

// Non-standard vec3-to-vec3 hash function.
vec3 hash33(vec3 p)
{
    float n = sin(dot(p, vec3(7, 157, 113)));    
    return fract(vec3(2097152, 262144, 32768)*n); 
}

// 2x2 matrix rotation.
mat2 rot2(float a)
{
	float c = cos(a); float s = sin(a);
	return mat2(c, s, -s, c);
}

// Tri-Planar blending function. Based on an old Nvidia tutorial.
vec3 tex3D( sampler2D tex, in vec3 p, in vec3 n )
{
    n = max((abs(n) - 0.2)*7., 0.001); // n = max(abs(n), 0.001), etc.
    n /= (n.x + n.y + n.z );   
	return (texture(tex, p.yz)*n.x + texture(tex, p.zx)*n.y + texture(tex, p.xy)*n.z).xyz;
}

// The triangle function that Shadertoy user Nimitz has used in various triangle noise demonstrations.
// See Xyptonjtroz - Very cool. Anyway, it's not really being used to its full potential here.
vec3 tri(in vec3 x){return abs(x-floor(x)-.5);} // Triangle function.

// The function used to perturb the walls of the cavern: There are infinite possibities, but this one is 
// just a cheap...ish routine - based on the triangle function - to give a subtle jaggedness. Not very fancy, 
// but it does a surprizingly good job at laying the foundations for a sharpish rock face. Obviously, more 
// layers would be more convincing. However, this is a GPU-draining distance function, so the finer details 
// are bump mapped.
float surfFunc(in vec3 p)
{
    // all have range: [0, 1]

#if defined(ENABLE_SURF_JAG1)
	return dot(tri(p*0.5 + tri(p*0.25 + 0.25).yzx), vec3(0.666));

#elif defined(ENABLE_SURF_JAG2)
    float n = dot(tri(p*0.5 + tri(p*0.25 + 0.25).yzx), vec3(0.444));
    p.xz = vec2(p.x + p.z, p.z - p.x) * 0.7071;
    return dot(tri(p*0.75 + tri(p*0.375 + 0.125).yzx), vec3(0.222)) + n; // Range [0, 1]

#elif defined(ENABLE_SURF_JAG3)
    return dot(tri(p*0.5 + tri(p*0.25).yzx), vec3(0.333)) + 
           sin(p.x*1.5+sin(p.y*2.+sin(p.z*2.5)))*0.25+0.25;

#elif defined(ENABLE_SURF_JAG4)
    return dot(tri(p*0.6 + tri(p*0.3).yzx), vec3(0.333)) + 
           sin(p.x*1.75+sin(p.y*2.+sin(p.z*2.25)))*0.25+0.25; // Range [0, 1]

#elif defined(ENABLE_SURF_JAG5)
    p *= 0.5;
    float n = dot(tri(p + tri(p*0.5).yzx), vec3(0.666*0.66));
    p *= 1.5;
    p.xz = vec2(p.x + p.z, p.z - p.x) * 1.7321*0.5;
    n += dot(tri(p + tri(p*0.5).yzx), vec3(0.666*0.34));
    return n;

#elif defined(ENABLE_SURF_JAG6)
    p *= 2.;
    float n = sin(p.x+sin(p.y+sin(p.z)))*0.57;
    p *= 1.5773;
    p.xy = vec2(p.x + p.y, p.y - p.x) * 1.7321*0.5;
    n += sin(p.x+sin(p.y+sin(p.z)))*0.28;
    p *= 1.5773;
    p.xy = vec2(p.x + p.y, p.y - p.x) * 1.7321*0.5;
    n += sin(p.x+sin(p.y+sin(p.z)))*0.15;
    return n*0.4+0.6;
#endif
}

// Cheap...ish smooth minimum function.
float smoothMinP( float a, float b, float smoothing )
{
    float h = clamp((b-a)*0.5/smoothing + 0.5, 0.0, 1.0 );
    return mix(b, a, h) - smoothing*h*(1.0-h);
}

// The path is a 2D sinusoid that varies over time, depending upon the frequencies, and amplitudes.
vec2 path(in float z){ float s = sin(z/24.)*cos(z/12.); return vec2(s*12., 0.); }

// Standard tunnel distance function with a bit of perturbation thrown into the mix. A winding 
// tunnel is just a tube with a smoothly shifting center as you traverse lengthwise. The walls 
// of the tunnels should be perturbed by some kind of 3D surface function... preferably a cheap 
// one with decent visual impact.
float map(vec3 p)
{
#if defined(ENABLE_ROUND_TUNNEL)
    // Round tunnel with floor using Euclidean distance: length(tun.xy)
    vec2 tun = abs(p.xy - path(p.z))*vec2(0.5, 0.7071);
    float n = 1.- length(tun.xy) + (0.5-surfFunc(p)); //max(tun.x, tun.y)
    return min(n, p.y + GROUND_HEIGHT);

#elif defined(ENABLE_SQUARE_TUNNEL)
    // Square tunnel using Chebyshev distance: max(abs(tun.x), abs(tun.y))
    vec2 tun = abs(p.xy - path(p.z))*vec2(0.5, 0.7071);
    //tun *= tun;
    float n = 1.- max(tun.x, tun.y) + (0.5-surfFunc(p));
    return min(n, p.y + GROUND_HEIGHT);

#elif defined(ENABLE_ROUNDED_SQUARE_TUNNEL)
    // Rounded square tunnel using Minkowski distance: pow(pow(abs(tun.x), n), pow(abs(tun.y), n), 1/n)
    vec2 tun = abs(p.xy - path(p.z))*vec2(0.5, 0.7071);
    tun = pow(tun, vec2(4.));
    float n =1.-pow(tun.x + tun.y, 1.0/4.) + (0.5-surfFunc(p));
    return min(n, p.y+GROUND_HEIGHT);

#elif defined(ENABLE_MINIMALISTS)
    // For the minimalists. :)
    float n = 0.5-surfFunc(p + 0.25*sign(p.y));
    n = min(GROUND_HEIGHT, 1. - sign(p.y)*n*0.75);
    return min(-p.y+n, p.y + n);
#endif
}

// Texture bump mapping. Four tri-planar lookups, or 12 texture lookups in total.
vec3 doBumpMap( sampler2D tex, in vec3 p, in vec3 nor, float bumpfactor )
{
    const float eps = 0.001;
    float ref = getGrey(tex3D(tex,  p , nor));                 
    vec3 grad = vec3( getGrey(tex3D(tex, vec3(p.x+eps, p.y, p.z), nor))-ref,
                      getGrey(tex3D(tex, vec3(p.x, p.y+eps, p.z), nor))-ref,
                      getGrey(tex3D(tex, vec3(p.x, p.y, p.z+eps), nor))-ref )/eps;
    grad -= nor*dot(nor, grad);
    return normalize( nor + grad*bumpfactor );
}

// Surface normal.
vec3 getNormal(in vec3 p)
{
	const float eps = 0.001;
	return normalize(vec3(
		map(vec3(p.x+eps,p.y,p.z))-map(vec3(p.x-eps,p.y,p.z)),
		map(vec3(p.x,p.y+eps,p.z))-map(vec3(p.x,p.y-eps,p.z)),
		map(vec3(p.x,p.y,p.z+eps))-map(vec3(p.x,p.y,p.z-eps))
	));
}

// Based on original by IQ.
float calculateAO(vec3 p, vec3 n)
{
    const float AO_SAMPLES = 5.0;
    float r = 0.0, w = 1.0, d;
    for (float i=1.0; i<AO_SAMPLES+1.1; i++){
        d = i/AO_SAMPLES;
        r += w*(d - map(p + n*d));
        w *= 0.5;
    }
    return 1.0-clamp(r,0.0,1.0);
}

// Cool curve function, by Shadertoy user, Nimitz.
// I think it's based on a discrete finite difference approximation to the continuous
// Laplace differential operator? Either way, it gives you the curvature of a surface, 
// which is pretty handy. I used it to do a bit of fake shadowing.
float curve(in vec3 p, in float w)
{
    vec2 e = vec2(-1., 1.)*w;
    float t1 = map(p + e.yxx), t2 = map(p + e.xxy);
    float t3 = map(p + e.xyx), t4 = map(p + e.yyy);
    return 0.125/(w*w) *(t1 + t2 + t3 + t4 - 4.*map(p));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
#ifdef ENABLE_AUTO_CAMERA
    //===================
    // original codes
    //===================

	// screen coordinates
	vec2 uv = (fragCoord - resolution.xy*0.5)/resolution.y;
	
	// camera setup
    vec3 lookAt = vec3(0.0, 0.0, time * CAMERA_SPEED); // "lookAt" position.
	vec3 camPos = lookAt + vec3(0.0, 0.1, -0.5); // Camera position, doubling as the ray origin.
 
    // light positioning. One is a little behind the camera, and the other is further down the tunnel.
 	vec3 light_pos = camPos + vec3(0.0, 0.125, -0.125); // Put it a bit in front of the camera.
	vec3 light_pos2 = camPos + vec3(0.0, 0.0, 6.0);     // Put it a bit in front of the camera.

	// Using the Z-value to perturb the XY-plane.
	// sending the "camera", "lookAt," and two lights down the tunnel.
    // The "path" function is synchronized with the distance function.
    // change to "path2" to traverse the other tunnel.
	lookAt.xy += path(lookAt.z);
	camPos.xy += path(camPos.z);
	light_pos.xy += path(light_pos.z);
	light_pos2.xy += path(light_pos2.z);

    // ray direction
    float FOV = PI/3.; // FOV - Field of view.
    vec3 forward = normalize(lookAt-camPos);
    vec3 right = normalize(vec3(forward.z, 0., -forward.x )); 
    vec3 up = cross(forward, right);
    vec3 rd = normalize(forward + FOV*uv.x*right + FOV*uv.y*up);

#else
    //===================
    // my codes... (to control view)
    //===================

    @import ./ray;
	
	// screen coordinates
    vec2 uv = screenPos;
	
	// camera setup
	vec3 lookAt = vec3(0.0, 0.0, time * CAMERA_SPEED); // "lookAt" position.
    vec3 camPos = -lookAt + rd; // Camera position, doubling as the ray origin.
 
    // light positioning. One is a little behind the camera, and the other is further down the tunnel.
 	vec3 light_pos = camPos - vec3(0.0, 0.125, -0.125);// Put it a bit in front of the camera.
	vec3 light_pos2 = camPos - vec3(0.0, 0.0, 6.0);// Put it a bit in front of the camera.

	// Using the Z-value to perturb the XY-plane.
	// sending the "camera", "lookAt," and two lights down the tunnel.
    // The "path" function is synchronized with the distance function.
    // change to "path2" to traverse the other tunnel.
	lookAt.xy += path(lookAt.z);
	camPos.xy += path(camPos.z);
	light_pos.xy += path(light_pos.z);
	light_pos2.xy += path(light_pos2.z);

#endif

	rd.xy *= rot2( -path(lookAt.z).x/32. );
    //rd.xz *= rot2( path(lookAt.z).x/32. );
		
    // Standard ray marching routine. I find that some system setups don't like anything other than
    // a "break" statement (by itself) to exit. 
	float t = 0.0, dt;
	for(int i=0; i<128; i++)
    {
		dt = map( camPos + rd * t );
		if( dt < 0.005 || t > 150.0 ){ break; } 
		t += dt * 0.75;
	}

	vec3 sceneCol = vec3(0.0);

	// The ray has effectively hit the surface, so light it up.
	if( dt < 0.005 )
    {
	    // The ray marching loop (above) exits when "dt" is less than a certain threshold, which in this 
        // case, is hardcoded to "0.005." However, the distance is still "dt" from the surface? By my logic, 
	    // adding the extra "dt" after breaking would gain a little more accuracy and effectively reduce 
	    // surface popping? Would that be correct? I tend to do this, but could be completely wrong, so if 
	    // someone could set me straight, it'd be appreciated. 
	    t += dt;

    	// Surface position and surface normal.
	    vec3 sp = t * rd+camPos;
	    vec3 sn = getNormal(sp);

        gl_FragDepth = getFragDepth( sp );

        // Texture scale factor.
        const float tSize0 = 1./2.; 
        const float tSize1 = 1./3.; 

    	// Texture-based bump mapping. Comment this line out to spoil the illusion.
	    if( sp.y < -(GROUND_HEIGHT - 0.005) )
            sn = doBumpMap( texChannel0, sp*tSize1, sn, 0.0125 );
	    else
            sn = doBumpMap( texChannel0, sp*tSize0, sn, 0.025 );

	    // Ambient occlusion.
	    float ao = calculateAO(sp, sn);

    	// Light direction vectors.
	    vec3 ld = light_pos-sp;
	    vec3 ld2 = light_pos2-sp;

        // Distance from respective lights to the surface point.
	    float distlpsp = max(length(ld), 0.001);
	    float distlpsp2 = max(length(ld2), 0.001);

    	// Normalize the light direction vectors.
	    ld /= distlpsp;
	    ld2 /= distlpsp2;

	    // Light attenuation, based on the distances above.
	    float atten = min(1./(distlpsp) + 1./(distlpsp2), 1.);

    	// Ambient light.
	    float ambience = 0.25;

    	// Diffuse lighting.
	    float diff = max( dot(sn, ld), 0.0);
	    float diff2 = max( dot(sn, ld2), 0.0);

    	// Specular lighting.
	    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.0 ), 8.);
	    float spec2 = pow(max( dot( reflect(-ld2, sn), -rd ), 0.0 ), 8.);

    	// Curvature.
	    float crv = clamp(curve(sp, 0.125)*0.5+0.5, .0, 1.);

	    // Fresnel term. Good for giving a surface a bit of a reflective glow.
        float fre = pow( clamp(dot(sn, rd) + 1., .0, 1.), 1.);

        // Obtaining the texel color. 
        vec3 texCol;
        if( sp.y < -(GROUND_HEIGHT - 0.005) )
        {
            texCol = tex3D( texChannel0, sp*tSize1, sn );
            texCol = texCol*0.5 + getGrey(texCol)*0.25 + 0.25;
        }
	    else texCol = tex3D( texChannel0, sp*tSize0, sn ); // Sandstone.

        // Jitter, if needed.
        //vec3 aniso = (0.5-hash33(sp))*fre*0.35;
	    //texCol = clamp(texCol + aniso, 0., 1.);

    	// Darkening the crevices. Otherwise known as cheap, scientifically-incorrect shadowing.	
	    float shading = crv * 0.5 + 0.5; 

    	// Combing the above terms to produce the final color. It was based more on acheiving a
        // certain aesthetic than science.

#if defined(ENABLE_GLOW)
        sceneCol = getGrey(texCol)*( (diff+diff2)*0.75 + ambience*0.25 ) + (spec + spec2)*texCol*1.5 + fre*crv*texCol.zyx*2.0;
#elif defined(ENABLE_SHINY)
        sceneCol = texCol*((diff+diff2)*vec3(1.0, 0.95, 0.9) + ambience + fre*fre*texCol) + (spec+spec2);
#endif
        //if( sp.y > -0.995 ) sceneCol += (spec + spec2)*texCol*0.5;

        // Shading
        sceneCol *= atten * shading * ao;

#ifdef ENABLE_EDGE_EFFECT
        //sceneCol *= clamp(1.-abs(curve(sp, 0.01)), .0, 1.); // add dense edges
        sceneCol *= pow(clamp(1.-(curve(sp, 0.0125)), 0., 1.), 0.5); // add thick edges
#endif
	}
    else gl_FragDepth = 0.99;

	fragColor = vec4( clamp(sceneCol, 0.0, 1.0), 1.0 );
}

void main()
{
    mainImage( gl_FragColor, gl_FragCoord.xy );
}