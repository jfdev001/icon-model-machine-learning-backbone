```{eval-rst}
:orphan:
```

(ref_tools_comin)=
# Community Interface (ComIn)

```{image} ComIn_logo_black.svg
:class: only-light
:align: right
:height: 100
:width: 200
```

```{image} ComIn_logo_white.svg
:class: only-dark
:align: right
:height: 100
:width: 200
```

Code contributions from different researchers and institutions ("third-party code") are usually not included in the main ICON code, but remain confined to project branches.
In any case, they add specific switches and calls to ICON's main loop, making the model code less readable.
On the other hand, to keep the third-party code compatible with new versions of ICON additional maintenance is required.
These problems are solved by providing a unified plugin interface: The ICON Community Interface ComIn.
While the core model remains unchanged, third-party code can be run alongside ICON, even if it is implemented in a programming language other than Fortran.
ComIn provides interfaces for **Fortran**, **C**, and **Python**.

The **Community Interface (ComIn)** organizes the data exchange and simulation events between the ICON model and "3rd party modules".
The concept can be logically divided into an Adapter Library and a Callback Register.

* **Adapter Library**: It is included in both, the ICON model and the 3rd party module.
It contains descriptive data structures, and regulates the access to existing and the creation of additional model variables.
* **Callback Register**: Subroutines of the 3rd party module may be called at pre-defined events during the model simulation.

::::{grid} 1 2 2 2
:gutter: 1 1 1 2

:::{grid-item-card}
:link: https://icon-comin.gitlab-pages.dkrz.de/comin/
**Documentation**
^^^
* User Guide
* Examples
:::

:::{grid-item-card}
:link: https://gitlab.dkrz.de/icon-comin/comin-training-exercises
**Exercises**
^^^
* Programming a Rather Simple ComIn Python Plugin
* Masking (Non-)Prognostic Cells
* Implementing a Diagnostic Function as a ComIn Plugin
:::

::::
