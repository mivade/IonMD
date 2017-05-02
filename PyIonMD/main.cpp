#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <ionmd/params.hpp>
#include <ionmd/simulation.hpp>
#include <ionmd/trap.hpp>

namespace py = pybind11;

using ionmd::SimParams;
using ionmd::Simulation;
using ionmd::SimStatus;
using ionmd::Trap;


PYBIND11_PLUGIN(ionmd)
{
    py::module m("ionmd", "IonMD Python bindings");

    py::enum_<SimStatus>(m, "Status")
        .value("IDLE", SimStatus::IDLE)
        .value("RUNNING", SimStatus::RUNNING)
        .value("FINISHED", SimStatus::FINISHED)
        .value("ERRORED", SimStatus::ERRORED);

    py::class_<SimParams>(m, "Params")
        .def(py::init())
        .def_readwrite("dt", &SimParams::dt)
        .def_readwrite("num_steps", &SimParams::num_steps)
        .def_readwrite("verbosity", &SimParams::verbosity)
        .def_readwrite("micromotion_enabled", &SimParams::micromotion_enabled)
        .def_readwrite("stochastic_enabled", &SimParams::stochastic_enabled)
        .def_readwrite("doppler_enabled", &SimParams::doppler_enabled)
        .def_readwrite("filename", &SimParams::filename)
        .def_readwrite("buffer_size", &SimParams::buffer_size);

    py::class_<Trap>(m, "Trap")
        .def(py::init())
        .def_readwrite("r0", &Trap::r0)
        .def_readwrite("z0", &Trap::z0)
        .def_readwrite("kappa", &Trap::kappa)
        .def_readwrite("omega_rf", &Trap::omega_rf)
        .def_readwrite("V_rf", &Trap::V_rf)
        .def_readwrite("U_dc", &Trap::U_dc)
        .def_readwrite("U_ec", &Trap::U_ec);

    py::class_<Simulation>(m, "Simulation")
        .def(py::init())
        .def_property("params", &Simulation::get_params, &Simulation::set_params)
        .def_property("trap", &Simulation::get_trap, &Simulation::set_trap)
        .def_readonly("status", &Simulation::status);
        .def("set_params", &Simulation::set_params)
        .def("set_trap", &Simulation::set_trap)
        .def("add_ion", &Simulation::add_ion)
        .def("start", &Simulation::start)

    return m.ptr();
}
