#ifndef IG_UGCORE_CONFIG_HPP
#define IG_UGCORE_CONFIG_HPP

#define FEATURE_REGISTRY_CLASS_NAME_MAP 0
#define FEATURE_REGISTRY_CLASS_GROUP_MAP 0

#define FEATURE_MEMORY_ALIGNED 0


#define FEATURE_MOVE_SYMMETRIC_MATRIX 0
/*
#if defined(FEATURE_MEMORY_ALIGNED) && FEATURE_MEMORY_ALIGNED == 1
alignas(64) value_type m_data[N][M];
#else
value_type m_data[N][M];
#endif
*/

#endif
