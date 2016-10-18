LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)
LOCAL_SRC_FILES := geometry.c om.c random.c algebra.c
LOCAL_C_INCLUDES := $(LOCAL_PATH)/include
LOCAL_MODULE     := openmotion
LOCAL_CFLAGS += -std=c99
include $(BUILD_SHARED_LIBRARY)


