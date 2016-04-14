#############################################################################
# Makefile for building: DLCA
# Generated by qmake (3.0) (Qt 5.5.1)
# Project:  ../DLCA/DLCA.pro
# Template: app
# Command: /Users/berntsim/Qt/5.5/clang_64/bin/qmake -spec macx-clang CONFIG+=x86_64 -o Makefile ../DLCA/DLCA.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang
CXX           = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++
DEFINES       = 
CFLAGS        = -pipe -O2 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.7 -Wall -W -fPIC $(DEFINES)
CXXFLAGS      = -pipe -O2 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -std=c++11 -stdlib=libc++ -mmacosx-version-min=10.7 -Wall -W -fPIC $(DEFINES)
INCPATH       = -I../DLCA -I. -I/usr/local/include -I../../../Qt/5.5/clang_64/mkspecs/macx-clang
QMAKE         = /Users/berntsim/Qt/5.5/clang_64/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = DLCA1.0.0
DISTDIR = /Users/berntsim/Documents/Master/build-DLCA-Desktop_Qt_5_5_1_clang_64bit-Release/.tmp/DLCA1.0.0
LINK          = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++
LFLAGS        = -headerpad_max_install_names -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -stdlib=libc++ -mmacosx-version-min=10.7
LIBS          = $(SUBLIBS) -L/usr/local/lib -lsfml-audio -lsfml-graphics -lsfml-system -lsfml-network -lsfml-window 
AR            = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ar cq
RANLIB        = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib -s
SED           = sed
STRIP         = 

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = ../DLCA/main.cpp \
		../DLCA/routines.cpp \
		../DLCA/containers.cpp 
OBJECTS       = main.o \
		routines.o \
		containers.o
DIST          = ../../../Qt/5.5/clang_64/mkspecs/features/spec_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/qdevice.pri \
		../../../Qt/5.5/clang_64/mkspecs/features/device_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/common/unix.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/macx.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/sanitize.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/gcc-base.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/gcc-base-mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/clang.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/clang-mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/qconfig.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bootstrap_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_clucene_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_platformsupport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webview_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/features/qt_functions.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/qt_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/macx-clang/qmake.conf \
		../../../Qt/5.5/clang_64/mkspecs/features/spec_post.prf \
		../DLCA/.qmake.stash \
		../../../Qt/5.5/clang_64/mkspecs/features/exclusive_builds.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/default_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/default_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/resolve_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/default_post.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/c++11.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/sdk.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/default_post.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/objective_c.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/warn_on.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/rez.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/testcase_targets.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/exceptions.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/yacc.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/lex.prf \
		../DLCA/DLCA.pro routines.h \
		containers.h ../DLCA/main.cpp \
		../DLCA/routines.cpp \
		../DLCA/containers.cpp
QMAKE_TARGET  = DLCA
DESTDIR       = #avoid trailing-slash linebreak
TARGET        = DLCA

####### Custom Compiler Variables
QMAKE_COMP_QMAKE_OBJECTIVE_CFLAGS = -pipe \
		-O2 \
		-isysroot \
		/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk \
		-std=c++11 \
		-stdlib=libc++ \
		-mmacosx-version-min=10.7 \
		-Wall \
		-W



first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: ../DLCA/DLCA.pro ../../../Qt/5.5/clang_64/mkspecs/macx-clang/qmake.conf ../../../Qt/5.5/clang_64/mkspecs/features/spec_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/qdevice.pri \
		../../../Qt/5.5/clang_64/mkspecs/features/device_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/common/unix.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/macx.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/sanitize.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/gcc-base.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/gcc-base-mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/clang.conf \
		../../../Qt/5.5/clang_64/mkspecs/common/clang-mac.conf \
		../../../Qt/5.5/clang_64/mkspecs/qconfig.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bootstrap_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_clucene_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_platformsupport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webview_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns.pri \
		../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri \
		../../../Qt/5.5/clang_64/mkspecs/features/qt_functions.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/qt_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/macx-clang/qmake.conf \
		../../../Qt/5.5/clang_64/mkspecs/features/spec_post.prf \
		.qmake.stash \
		../../../Qt/5.5/clang_64/mkspecs/features/exclusive_builds.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/default_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/default_pre.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/resolve_config.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/default_post.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/c++11.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/sdk.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/default_post.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/objective_c.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/warn_on.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/mac/rez.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/testcase_targets.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/exceptions.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/yacc.prf \
		../../../Qt/5.5/clang_64/mkspecs/features/lex.prf \
		../DLCA/DLCA.pro
	$(QMAKE) -spec macx-clang CONFIG+=x86_64 -o Makefile ../DLCA/DLCA.pro
../../../Qt/5.5/clang_64/mkspecs/features/spec_pre.prf:
../../../Qt/5.5/clang_64/mkspecs/qdevice.pri:
../../../Qt/5.5/clang_64/mkspecs/features/device_config.prf:
../../../Qt/5.5/clang_64/mkspecs/common/unix.conf:
../../../Qt/5.5/clang_64/mkspecs/common/mac.conf:
../../../Qt/5.5/clang_64/mkspecs/common/macx.conf:
../../../Qt/5.5/clang_64/mkspecs/common/sanitize.conf:
../../../Qt/5.5/clang_64/mkspecs/common/gcc-base.conf:
../../../Qt/5.5/clang_64/mkspecs/common/gcc-base-mac.conf:
../../../Qt/5.5/clang_64/mkspecs/common/clang.conf:
../../../Qt/5.5/clang_64/mkspecs/common/clang-mac.conf:
../../../Qt/5.5/clang_64/mkspecs/qconfig.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcollision_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dcore_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dinput_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dlogic_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquick_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3dquickrenderer_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_3drenderer_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bluetooth_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_bootstrap_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_clucene_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_concurrent_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_core_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_dbus_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_declarative_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designer_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_designercomponents_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_enginio_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_gui_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_help_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_location_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_macextras_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimedia_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_multimediawidgets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_network_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_nfc_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_opengl_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_openglextensions_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_platformsupport_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_positioning_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_printsupport_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qml_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmldevtools_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qmltest_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_qtmultimediaquicktools_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quick_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickparticles_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_quickwidgets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_script_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_scripttools_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sensors_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_serialport_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_sql_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_svg_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_testlib_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uiplugin.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_uitools_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webchannel_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webengine_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginecore_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webenginewidgets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkit_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webkitwidgets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_websockets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_webview_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_widgets_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xml_private.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns.pri:
../../../Qt/5.5/clang_64/mkspecs/modules/qt_lib_xmlpatterns_private.pri:
../../../Qt/5.5/clang_64/mkspecs/features/qt_functions.prf:
../../../Qt/5.5/clang_64/mkspecs/features/qt_config.prf:
../../../Qt/5.5/clang_64/mkspecs/macx-clang/qmake.conf:
../../../Qt/5.5/clang_64/mkspecs/features/spec_post.prf:
.qmake.stash:
../../../Qt/5.5/clang_64/mkspecs/features/exclusive_builds.prf:
../../../Qt/5.5/clang_64/mkspecs/features/default_pre.prf:
../../../Qt/5.5/clang_64/mkspecs/features/mac/default_pre.prf:
../../../Qt/5.5/clang_64/mkspecs/features/resolve_config.prf:
../../../Qt/5.5/clang_64/mkspecs/features/default_post.prf:
../../../Qt/5.5/clang_64/mkspecs/features/c++11.prf:
../../../Qt/5.5/clang_64/mkspecs/features/mac/sdk.prf:
../../../Qt/5.5/clang_64/mkspecs/features/mac/default_post.prf:
../../../Qt/5.5/clang_64/mkspecs/features/mac/objective_c.prf:
../../../Qt/5.5/clang_64/mkspecs/features/warn_on.prf:
../../../Qt/5.5/clang_64/mkspecs/features/mac/rez.prf:
../../../Qt/5.5/clang_64/mkspecs/features/testcase_targets.prf:
../../../Qt/5.5/clang_64/mkspecs/features/exceptions.prf:
../../../Qt/5.5/clang_64/mkspecs/features/yacc.prf:
../../../Qt/5.5/clang_64/mkspecs/features/lex.prf:
../DLCA/DLCA.pro:
qmake: FORCE
	@$(QMAKE) -spec macx-clang CONFIG+=x86_64 -o Makefile ../DLCA/DLCA.pro

qmake_all: FORCE


all: Makefile $(TARGET)

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


distclean: clean 
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) .qmake.stash
	-$(DEL_FILE) Makefile


####### Sub-libraries

check: first

compiler_objective_c_make_all:
compiler_objective_c_clean:
compiler_rez_source_make_all:
compiler_rez_source_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: ../DLCA/main.cpp /usr/local/include/SFML/Graphics.hpp \
		/usr/local/include/SFML/Window.hpp \
		/usr/local/include/SFML/System.hpp \
		/usr/local/include/SFML/Config.hpp \
		/usr/local/include/SFML/System/Clock.hpp \
		/usr/local/include/SFML/System/Export.hpp \
		/usr/local/include/SFML/System/Time.hpp \
		/usr/local/include/SFML/System/Err.hpp \
		/usr/local/include/SFML/System/FileInputStream.hpp \
		/usr/local/include/SFML/System/InputStream.hpp \
		/usr/local/include/SFML/System/NonCopyable.hpp \
		/usr/local/include/SFML/System/Lock.hpp \
		/usr/local/include/SFML/System/MemoryInputStream.hpp \
		/usr/local/include/SFML/System/Mutex.hpp \
		/usr/local/include/SFML/System/Sleep.hpp \
		/usr/local/include/SFML/System/String.hpp \
		/usr/local/include/SFML/System/Utf.hpp \
		/usr/local/include/SFML/System/Utf.inl \
		/usr/local/include/SFML/System/String.inl \
		/usr/local/include/SFML/System/Thread.hpp \
		/usr/local/include/SFML/System/Thread.inl \
		/usr/local/include/SFML/System/ThreadLocal.hpp \
		/usr/local/include/SFML/System/ThreadLocalPtr.hpp \
		/usr/local/include/SFML/System/ThreadLocalPtr.inl \
		/usr/local/include/SFML/System/Vector2.hpp \
		/usr/local/include/SFML/System/Vector2.inl \
		/usr/local/include/SFML/System/Vector3.hpp \
		/usr/local/include/SFML/System/Vector3.inl \
		/usr/local/include/SFML/Window/Context.hpp \
		/usr/local/include/SFML/Window/Export.hpp \
		/usr/local/include/SFML/Window/GlResource.hpp \
		/usr/local/include/SFML/Window/ContextSettings.hpp \
		/usr/local/include/SFML/Window/Event.hpp \
		/usr/local/include/SFML/Window/Joystick.hpp \
		/usr/local/include/SFML/Window/Keyboard.hpp \
		/usr/local/include/SFML/Window/Mouse.hpp \
		/usr/local/include/SFML/Window/Sensor.hpp \
		/usr/local/include/SFML/Window/Touch.hpp \
		/usr/local/include/SFML/Window/VideoMode.hpp \
		/usr/local/include/SFML/Window/Window.hpp \
		/usr/local/include/SFML/Window/WindowHandle.hpp \
		/usr/local/include/SFML/Window/WindowStyle.hpp \
		/usr/local/include/SFML/Graphics/BlendMode.hpp \
		/usr/local/include/SFML/Graphics/Export.hpp \
		/usr/local/include/SFML/Graphics/CircleShape.hpp \
		/usr/local/include/SFML/Graphics/Shape.hpp \
		/usr/local/include/SFML/Graphics/Drawable.hpp \
		/usr/local/include/SFML/Graphics/RenderStates.hpp \
		/usr/local/include/SFML/Graphics/Transform.hpp \
		/usr/local/include/SFML/Graphics/Rect.hpp \
		/usr/local/include/SFML/Graphics/Rect.inl \
		/usr/local/include/SFML/Graphics/Transformable.hpp \
		/usr/local/include/SFML/Graphics/VertexArray.hpp \
		/usr/local/include/SFML/Graphics/Vertex.hpp \
		/usr/local/include/SFML/Graphics/Color.hpp \
		/usr/local/include/SFML/Graphics/PrimitiveType.hpp \
		/usr/local/include/SFML/Graphics/ConvexShape.hpp \
		/usr/local/include/SFML/Graphics/Font.hpp \
		/usr/local/include/SFML/Graphics/Glyph.hpp \
		/usr/local/include/SFML/Graphics/Texture.hpp \
		/usr/local/include/SFML/Graphics/Image.hpp \
		/usr/local/include/SFML/Graphics/RectangleShape.hpp \
		/usr/local/include/SFML/Graphics/RenderTarget.hpp \
		/usr/local/include/SFML/Graphics/View.hpp \
		/usr/local/include/SFML/Graphics/RenderTexture.hpp \
		/usr/local/include/SFML/Graphics/RenderWindow.hpp \
		/usr/local/include/SFML/Graphics/Shader.hpp \
		/usr/local/include/SFML/Graphics/Sprite.hpp \
		/usr/local/include/SFML/Graphics/Text.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o ../DLCA/main.cpp

routines.o: ../DLCA/routines.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o routines.o ../DLCA/routines.cpp

containers.o: ../DLCA/containers.cpp ../DLCA/containers.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o containers.o ../DLCA/containers.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:
