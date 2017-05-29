/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_keys
#define __H__UG_keys

namespace ug{
namespace promesh{

enum ModifierKeys{
	MK_NONE		= 0,
	MK_CTRL		= 1,
	MK_SHIFT	= 2,
	MK_ALT		= 4
};

enum ShortcutKeys {
	Key_Escape	= 0x01000000,
	Key_Tab	= 0x01000001,
	Key_Backtab	= 0x01000002,
	Key_Backspace	= 0x01000003,
	Key_Return	= 0x01000004,
	Key_Enter	= 0x01000005,
	Key_Insert	= 0x01000006,
	Key_Delete	= 0x01000007,
	Key_Pause	= 0x01000008,
	Key_Print	= 0x01000009,
	Key_SysReq	= 0x0100000a,
	Key_Clear	= 0x0100000b,
	Key_Home	= 0x01000010,
	Key_End	= 0x01000011,
	Key_Left	= 0x01000012,
	Key_Up	= 0x01000013,
	Key_Right	= 0x01000014,
	Key_Down	= 0x01000015,
	Key_PageUp	= 0x01000016,
	Key_PageDown	= 0x01000017,
	Key_Shift	= 0x01000020,
	Key_Control	= 0x01000021,
	Key_Meta	= 0x01000022,
	Key_Alt	= 0x01000023,
	Key_AltGr	= 0x01001103,
	Key_CapsLock	= 0x01000024,
	Key_NumLock	= 0x01000025,
	Key_ScrollLock	= 0x01000026,
	Key_F1	= 0x01000030,
	Key_F2	= 0x01000031,
	Key_F3	= 0x01000032,
	Key_F4	= 0x01000033,
	Key_F5	= 0x01000034,
	Key_F6	= 0x01000035,
	Key_F7	= 0x01000036,
	Key_F8	= 0x01000037,
	Key_F9	= 0x01000038,
	Key_F10	= 0x01000039,
	Key_F11	= 0x0100003a,
	Key_F12	= 0x0100003b,
	Key_F13	= 0x0100003c,
	Key_F14	= 0x0100003d,
	Key_F15	= 0x0100003e,
	Key_F16	= 0x0100003f,
	Key_F17	= 0x01000040,
	Key_F18	= 0x01000041,
	Key_F19	= 0x01000042,
	Key_F20	= 0x01000043,
	Key_F21	= 0x01000044,
	Key_F22	= 0x01000045,
	Key_F23	= 0x01000046,
	Key_F24	= 0x01000047,
	Key_F25	= 0x01000048,
	Key_F26	= 0x01000049,
	Key_F27	= 0x0100004a,
	Key_F28	= 0x0100004b,
	Key_F29	= 0x0100004c,
	Key_F30	= 0x0100004d,
	Key_F31	= 0x0100004e,
	Key_F32	= 0x0100004f,
	Key_F33	= 0x01000050,
	Key_F34	= 0x01000051,
	Key_F35	= 0x01000052,
	Key_Super_L	= 0x01000053,
	Key_Super_R	= 0x01000054,
	Key_Menu	= 0x01000055,
	Key_Hyper_L	= 0x01000056,
	Key_Hyper_R	= 0x01000057,
	Key_Help	= 0x01000058,
	Key_Direction_L	= 0x01000059,
	Key_Direction_R	= 0x01000060,
	Key_Space	= 0x20,
	Key_Any		= 0x20,
	Key_Exclam	= 0x21,
	Key_QuoteDbl	= 0x22,
	Key_NumberSign	= 0x23,
	Key_Dollar	= 0x24,
	Key_Percent	= 0x25,
	Key_Ampersand	= 0x26,
	Key_Apostrophe	= 0x27,
	Key_ParenLeft	= 0x28,
	Key_ParenRight	= 0x29,
	Key_Asterisk	= 0x2a,
	Key_Plus	= 0x2b,
	Key_Comma	= 0x2c,
	Key_Minus	= 0x2d,
	Key_Period	= 0x2e,
	Key_Slash	= 0x2f,
	Key_0	= 0x30,
	Key_1	= 0x31,
	Key_2	= 0x32,
	Key_3	= 0x33,
	Key_4	= 0x34,
	Key_5	= 0x35,
	Key_6	= 0x36,
	Key_7	= 0x37,
	Key_8	= 0x38,
	Key_9	= 0x39,
	Key_Colon	= 0x3a,
	Key_Semicolon	= 0x3b,
	Key_Less	= 0x3c,
	Key_Equal	= 0x3d,
	Key_Greater	= 0x3e,
	Key_Question	= 0x3f,
	Key_At	= 0x40,
	Key_A	= 0x41,
	Key_B	= 0x42,
	Key_C	= 0x43,
	Key_D	= 0x44,
	Key_E	= 0x45,
	Key_F	= 0x46,
	Key_G	= 0x47,
	Key_H	= 0x48,
	Key_I	= 0x49,
	Key_J	= 0x4a,
	Key_K	= 0x4b,
	Key_L	= 0x4c,
	Key_M	= 0x4d,
	Key_N	= 0x4e,
	Key_O	= 0x4f,
	Key_P	= 0x50,
	Key_Q	= 0x51,
	Key_R	= 0x52,
	Key_S	= 0x53,
	Key_T	= 0x54,
	Key_U	= 0x55,
	Key_V	= 0x56,
	Key_W	= 0x57,
	Key_X	= 0x58,
	Key_Y	= 0x59,
	Key_Z	= 0x5a,
	Key_BracketLeft	= 0x5b,
	Key_Backslash	= 0x5c
};

}	
}//	end of namespace

#endif	//__H__UG_keys
