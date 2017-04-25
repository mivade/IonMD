#include <pybind11/pybind11.h>

#include <ionmd/params.hpp>
#include <ionmd/simulation.hpp>
#include <ionmd/ion.hpp>
#include <ionmd/trap.hpp>

namespace py = pybind11;

using ionmd::SimParams;
using ionmd::Simulation;
using ionmd::Ion;
using ionmd::Trap;


PYBIND11_PLUGIN(pyionmd)
{
    py::module m("pyionmd", "IonMD Python bindings");

    py::class_<SimParams>(m, "SimParams")
        .def(py::init())
        .def_readwrite("dt", &SimParams::dt)
        .def_readwrite("num_steps", &SimParams::num_steps)
        .def_readwrite("verbosity", &SimParams::verbosity)
        .def_readwrite("micromotion_enabled", &SimParams::micromotion_enabled)
        .def_readwrite("stochastic_enabled", &SimParams::stochastic_enabled)
        .def_readwrite("doppler_enabled", &SimParams::doppler_enabled)
        .def_readwrite("filename", &SimParams::filename)
        .def_readwrite("buffer_size", &SimParams::buffer_size);

    py::class_<Ion>(m, "Ion")
        .def(py::init<SimParams &p, Trap &t, const double m, const double Z>);

    // py::class_<Trap>(m, "Trap")
    //     .def(py::init());

    py::class_<Simulation>(m, "Simulation")
        .def(py::init())
        .def("set_params", &Simulation::set_params)
        .def("set_trap", &Simulation::set_trap)
//        .def("make_ion" &Simulation::make_ion)
        .def("set_ions", &Simulation::set_ions)
        .def("run", &Simulation::run);

    return m.ptr();
}
