uniform sampler2D t_pos;

void main()
{
  vec4 pos = texture( t_pos, position.xy );
  float dist = length(cameraPosition - pos.xyz);

  //gl_PointSize = clamp( 500.0 / dist, 10.0, 50.0 );
  //gl_PointSize = max( 50.0 * exp(-0.005 * dist), 10.0 );
  gl_PointSize = max( 20.0 - 0.1 * dist, 5.0 );

  gl_Position = projectionMatrix * modelViewMatrix * vec4( pos.xyz, 1.0 );
}