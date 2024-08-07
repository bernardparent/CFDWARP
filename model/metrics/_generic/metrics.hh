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

typedef struct {
  int METRICSMODEL;
  double rmin_axisymmetric;
} gl_model_metrics_t;
