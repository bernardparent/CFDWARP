#define _METRICS_GENERIC
#define _METRICS_METHOD "Generic"
#define _METRICS_ACTIONNAME "Metrics"

#define METRICSMODEL_VIVIANDVINOKUR 1
#ifdef _3D
  #define METRICSMODEL_FREESTREAMPRESERVING 2
#endif
#ifdef _2D
  #define METRICSMODEL_AXISYMMETRIC 3
#endif

#define METRICSRMIN_OFF 0
#define METRICSRMIN_ON 1

typedef struct {
  int METRICSMODEL;
#ifdef _2D
  double axisymmetric_min_radius,axisymmetric_slice_angle;
#endif
} gl_model_metrics_t;
