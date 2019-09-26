#ifndef RAYMARCH_CURL
#define RAYMARCH_CURL

vec3 snoiseVec3( vec3 x )
{
  float s  = snoise(vec3( x ));
  float s1 = snoise(vec3( x.y - 19.1 , x.z + 33.4 , x.x + 47.2 ));
  float s2 = snoise(vec3( x.z + 74.2 , x.x - 124.5 , x.y + 99.4 ));
  vec3 c = vec3( s , s1 , s2 );
  return c;
}

vec3 curlNoise( vec3 p )
{
  const float e = 1e-1;
  vec3 dx = vec3( e   , 0.0 , 0.0 );
  vec3 dy = vec3( 0.0 , e   , 0.0 );
  vec3 dz = vec3( 0.0 , 0.0 , e   );
  vec3 p_x0 = snoiseVec3( p - dx );
  vec3 p_x1 = snoiseVec3( p + dx );
  vec3 p_y0 = snoiseVec3( p - dy );
  vec3 p_y1 = snoiseVec3( p + dy );
  vec3 p_z0 = snoiseVec3( p - dz );
  vec3 p_z1 = snoiseVec3( p + dz );
  float x = p_y1.z - p_y0.z - p_z1.y + p_z0.y;
  float y = p_z1.x - p_z0.x - p_x1.z + p_x0.z;
  float z = p_x1.y - p_x0.y - p_y1.x + p_y0.x;
  const float divisor = 1.0 / ( 2.0 * e );
  return normalize( vec3( x , y , z ) * divisor );
}

vec4 snoiseVec4( vec4 x )
{
  float s  = snoise(vec4( x ));
  float s1 = snoise(vec4( x.y - 19.1 , x.z + 33.4 , x.w + 47.2 , x.x + 12.2 ));
  float s2 = snoise(vec4( x.z + 74.2 , x.w - 124.5 , x.x + 99.4 , x.y - 123.2 ));
  float s3 = snoise(vec4( x.w + 21.2 , x.x - 52.5 , x.y + 60.4 , x.z + 42.2 ));
  vec4 c = vec4( s , s1 , s2 , s3 );
  return c;
}

vec4 curlNoise( vec4 p )
{
  const float e = 0.1;
  vec4 dx = vec4( e   , 0.0 , 0.0 , 0.0 );
  vec4 dy = vec4( 0.0 , e   , 0.0 , 0.0 );
  vec4 dz = vec4( 0.0 , 0.0 , e   , 0.0 );
  vec4 dw = vec4( 0.0 , 0.0 , 0.0 , e  );
  vec4 p_x0 = snoiseVec4( p - dx );
  vec4 p_x1 = snoiseVec4( p + dx );
  vec4 p_y0 = snoiseVec4( p - dy );
  vec4 p_y1 = snoiseVec4( p + dy );
  vec4 p_z0 = snoiseVec4( p - dz );
  vec4 p_z1 = snoiseVec4( p + dz );
  vec4 p_w0 = snoiseVec4( p - dw );
  vec4 p_w1 = snoiseVec4( p + dw );
  float x = p_y1.z - p_y0.z - p_z1.w + p_z0.w + p_w0.y - p_w1.y;
  float y = p_z1.w - p_z0.w - p_w1.x + p_w0.x + p_x0.z - p_x1.z;
  float z = p_w1.x - p_w0.x - p_x1.y + p_x0.y + p_y0.w - p_y1.w;
  float w = p_x1.y - p_x0.y - p_y1.z + p_y0.z + p_z0.x - p_y1.x;
  const float divisor = 1.0 / ( 2.0 * e );
  return normalize( vec4( x , y , z, w ) * divisor );
}

#endif // RAYMARCH_CURL