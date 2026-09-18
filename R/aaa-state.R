# Internal mutable state (namespace bindings are locked; this env is not).
.epito <- new.env(parent = emptyenv())
.epito$netmhcpan_use_wsl <- FALSE
