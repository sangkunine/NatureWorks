uniform sampler2D t_oPos;
uniform sampler2D t_pos;
uniform vec2 resolution;

uniform float dTime;
uniform float curTime;
uniform vec3 emitPos;       // emitter position (x,y,z)
uniform vec3 emitVec;       // emit direction and its length (not needed for radialCurl)
uniform float curlType;     // 0=radialCurl, 1=directionalCurl
uniform float curlShape;    // shape of curl (0 ~ 1) ==> 0.0(rocketJet), 0.5(smoking), 1.0(floating)
uniform float curlSpeed;    // moving speed of particles (0 ~ 1)
uniform float curlSpread;   // level of spread of particles (0 ~ 1)

in vec2 vUv;

@import ./simplex;
@import ./curl;

float rand(vec2 co)
{
  return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

void radialCurl()
{
  //vec3 emitPos = vec3(0.0);

  vec2 uv = gl_FragCoord.xy / resolution;
  //vec4 oPos = texture( t_oPos, uv );
  vec4 pos  = texture( t_pos, uv );

  float shape = mix(0.04, 0.01, curlShape);
  vec3 curl = curlNoise(vec4(pos.xyz * shape, curTime * 1.1 * 2.5)).xyz;
  float speed = mix(0.5, 2.5, curlSpeed);
  vec3 vel = curl * speed;

  vec4 newPos;
  newPos.xyz = pos.xyz + vel;
  newPos.w = pos.w + dTime; // pos.w = particle's time of life (age)

  if( newPos.w > 10.0 ) // died at 10(secs) old
  {
    newPos.xyz = emitPos;
    float spread = mix(0.1, 200.0, curlSpread);
    newPos.x += spread*(rand(uv*curTime+0.7)-0.5);
    newPos.y += spread*(rand(uv*curTime+1.3)-0.5);
    newPos.z += spread*(rand(uv*curTime+2.7)-0.5);
    newPos.w = 10.0*rand(uv*curTime+3.9);
  }

  gl_FragColor = vec4( newPos.xyz, newPos.w );
}

void directionalCurl()
{
  float emitLength = length(emitVec);
  vec3 emitDir = emitVec / emitLength;

  //vec3 curlParams = vec3(0.5, 0.4, 0.0);// <== default
  float shape = 0.5 * exp( -curlShape * 4.6 );
  float speed = 0.5 + 5.0 * curlSpeed;
  float spread = 0.1 + 2.0 * curlSpread;

  vec2 uv = gl_FragCoord.xy / resolution;
  //vec4 oPos = texture( t_oPos, uv );
  vec4 pos  = texture( t_pos, uv );

  vec3 curl = curlNoise(vec4(pos.xyz * shape, curTime * 0.5 * speed)).xyz;
  vec3 vel = curl*0.1;

  float t = dot(emitDir, pos.xyz - emitPos);
  float v = pow(smoothstep(0.0, emitLength, t), 2.0);
  v = 0.05 + v * (0.2 + rand(uv) * 0.2);
  v = clamp(v, 0.01, 5.0);
  vel += emitDir * v;
  vel *= speed;
  vec3 newPos = pos.xyz + vel * dTime*400.0;

  t = dot(emitDir, newPos - emitPos);
  if( t > emitLength ) {
    newPos = emitPos;
    newPos.x += emitLength*spread*(rand(uv + vec2(21.3, 63.21))-0.5);
    newPos.y += emitLength*spread*(rand(uv + vec2(32.3, 734.21))-0.5);
    newPos.z += emitLength*spread*(rand(uv + vec2(127.3, 31.21))-0.5);

    t = dot(emitDir, newPos - emitPos);
    if( t < 0.0 ) { newPos = newPos - t*emitDir; t = 0.0; }
    vec3 tpos = emitPos + t*emitDir;
    vec3 rvec = newPos - tpos;
    float r = (1.0 - cos(smoothstep(0.0, 0.2*emitLength, t) * 3.141592654)) * 0.5;
    newPos = newPos + (r-1.0) * rvec;
  }

  gl_FragColor = vec4( newPos, 1.0 );
}

void main()
{
  if( curlType < 0.5 )
    radialCurl();       // curlType = 0.0
  else
    directionalCurl();  // curlType = 1.0
}