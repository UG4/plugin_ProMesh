/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_promesh_registry
#define __H__UG_promesh_registry

#include <map>
#include <set>
#include <string>
#include <vector>
#include "keys.h"
#include "registry/registry.h"
#include "boost/mpl/assert.hpp"
#include "boost/mpl/contains.hpp"
namespace ug{
namespace promesh{

/// \addtogroup promesh
/// \{
enum RegistryTargets{
	RT_NONE		= 0,
	RT_UGSCRIPT	= 1,
	RT_PROMESH	= 1 << 1,

	RT_DEFAULT = RT_UGSCRIPT | RT_PROMESH,
	RT_NO_PROMESH = RT_UGSCRIPT,
	RT_NO_UGSCRIPT = RT_PROMESH
};

namespace detail {

///	All functions registered in the ProMeshRegistry are encapsulated in a ProMeshFunction
/**	shortcut keys are specified as enumerated in ug::promesh::ShortcutKeys*/
class UG_API ProMeshFunction{
	public:
		ProMeshFunction (	bridge::ExportedFunction* exportedFunction,
							int priority,
							int groupPriority,
							unsigned int target,
							int shortcutKey = 0,
							ModifierKeys modifierKey = MK_NONE)
		:
			m_exportedFunction(exportedFunction),
			m_priority(priority),
			m_groupPriority(groupPriority),
			m_target(target),
			m_shortcutKey(shortcutKey),
			m_modifierKey(modifierKey)
		{}

		bridge::ExportedFunction*
		exported_function()								{return m_exportedFunction;}

		const bridge::ExportedFunction*
		exported_function() const						{return m_exportedFunction;}

		int priority() const							{return m_priority;}
		int group_priority() const						{return m_groupPriority;}

		bool has_target(RegistryTargets t) const 		{return (m_target & t) == t;}

		bool operator < (const ProMeshFunction& f) const
		{
			if(m_priority != f.priority())
				return m_priority < f.priority();
			if(m_groupPriority != f.group_priority())
				return m_groupPriority < f.group_priority();
			if(m_exportedFunction->group() != f.exported_function()->group())
				return m_exportedFunction->group() < f.exported_function()->group();
			if(m_exportedFunction->name() != f.exported_function()->name())
				return m_exportedFunction->name() < f.exported_function()->name();
			return m_exportedFunction->num_parameter() < f.exported_function()->num_parameter();
		}

		int shortcut_key() const					{return m_shortcutKey;}
		ModifierKeys shortcut_modifier_key() const	{return m_modifierKey;}

	private:
		bridge::ExportedFunction*	m_exportedFunction;
		int							m_priority;
		int							m_groupPriority;
		unsigned int				m_target;
		int							m_shortcutKey;
		ModifierKeys				m_modifierKey;
};

}//	end of namespace detail

///	Register functions for ug-script and ProMesh through this class
/** The ProMeshRegistry is a small wrapper for ug::bridge::Registry and allows
 * for the specification of additional parameters, which are later on used by
 * ProMesh to automatically generate tools for all registered functions.*/
class UG_API ProMeshRegistry{
	public:
		typedef std::multiset<detail::ProMeshFunction>		ProMeshFunctionSet;
		typedef ProMeshFunctionSet::iterator				func_iter_t;
		typedef ProMeshFunctionSet::const_iterator			const_func_iter_t;

		ProMeshRegistry(bridge::Registry* reg) : m_reg(reg), m_counter(0) {}

	/**
	 * @brief adds a function to the registry
	 * @param funcName the name of the function
	 * @param func function pointer of the function
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param retValInfos string documenting what the function returns (optional)
	 * @param paramInfos string documenting the parameters of the function
	 * seperate parameters with an # e.g. "x#y#z" (don't specify the type of the values)  (optional)
	 * @param toolTip small documentation for the function (optional)
	 * @param help help string for the function
	 * @sa \ref pageUG4Registry
	 * @sa \ref secSTHowToSpecifyParameterInformation
	 *
	 * References the template function proxy_function<TFunc> and stores
	 * it with the FuntionWrapper.
	 */
		template<typename TFunc>
		ProMeshRegistry& add_function(	std::string funcName,
										TFunc func,
										std::string group = "",
										std::string retValInfos = "",
										std::string paramInfos = "",
										std::string tooltip = "",
										std::string help = "",
										unsigned int target = RT_DEFAULT,
										int shortcutKey = 0,
										ModifierKeys modifyerKey = MK_NONE)
		{
			using namespace bridge;
			ExportedFunction* ef =
				m_reg->add_and_get_function(
							funcName, func, group, retValInfos,
							paramInfos, tooltip, help);
			
			int& groupPriority = m_groupPriority[group];
			if(!groupPriority){
				groupPriority = m_counter;
			}

			m_funcSet.insert(
					detail::ProMeshFunction(ef, m_counter, groupPriority,
											target, shortcutKey, modifyerKey));
			++m_counter;
			return *this;
		}

	/**
	 * @brief Register a class at this registry
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass>
		bridge::ExportedClass<TClass>&
		add_class_(	std::string className,
		            std::string group = "",
		            std::string tooltip = "")
		{
			return m_reg->add_class_<TClass>(className, group, tooltip);
		}

	/**
	 * @brief Register a class at this registry together with its base class
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass, typename TBaseClass>
		bridge::ExportedClass<TClass>&
		add_class_(	std::string className,
		            std::string group = "",
		            std::string tooltip = "")
		{
			return m_reg->add_class_<TClass, TBaseClass>(className, group, tooltip);
		}

	/**
	 * @brief Register a class at this registry together with two base classes
	 * @param className name of the class to appear in the registry
	 * @param group registry group. use / for subgroups e.g. ug4/mygroup/mysubgroup (optional)
	 * @param toolTip describing text for the class (optional)
	 */
		template <typename TClass, typename TBaseClass1, typename TBaseClass2>
		bridge::ExportedClass<TClass>&
		add_class_(	std::string className,
		            std::string group = "",
		            std::string tooltip = "")
		{
			return m_reg->add_class_<TClass, TBaseClass1, TBaseClass2>(
													className, group, tooltip);
		}

	///	returns a pointer to the ug::bridge::Registry which is encapsulated by this class.
		bridge::Registry* registry()				{return m_reg;}


		func_iter_t functions_begin()				{return m_funcSet.begin();}
		const_func_iter_t functions_begin() const	{return m_funcSet.begin();}
		func_iter_t functions_end()					{return m_funcSet.end();}
		const_func_iter_t functions_end() const		{return m_funcSet.end();}

	private:
		bridge::Registry* 			m_reg;
		ProMeshFunctionSet			m_funcSet;
		int							m_counter;
		std::map<std::string, int>	m_groupPriority;
};

/// \}

}//	end of namespace
}//	end of namespace

#endif	//__H__UG_promesh_registry
