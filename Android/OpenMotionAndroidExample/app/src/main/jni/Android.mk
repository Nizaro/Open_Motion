LOCAL_PATH := $(call my-dir)
OPENMOTION_SDK_PATH := /home/thomas/These/Library/OpenMotion/android-native

include $(CLEAR_VARS)

include $(OPENMOTION_SDK_PATH)/jni/OpenMotion.mk


LOCAL_CFLAGS += -I$(OPENMOTION_SDK_PATH)/jni/include

LOCAL_MODULE    := native_lib
LOCAL_SRC_FILES := main_native.c
LOCAL_C_INCLUDE := $(LOCAL_PATH)
LOCAL_LDLIBS    += -lm -llog -landroid -lz -lEGL -lGLESv2
LOCAL_CFLAGS += -std=c99
LOCAL_SHARED_LIBRARIES += openmotion

include $(BUILD_SHARED_LIBRARY)
