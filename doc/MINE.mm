<map version="1.0.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1393941195253" ID="ID_839841935" MODIFIED="1393941202711" TEXT="MINE">
<node CREATED="1393943512085" ID="ID_151805542" MODIFIED="1393943515808" POSITION="right" TEXT="Algorithms">
<node CREATED="1393943515810" ID="ID_265344647" MODIFIED="1393943863418" TEXT="MaxMI">
<font BOLD="true" NAME="SansSerif" SIZE="12"/>
<node CREATED="1393943522666" ID="ID_338849859" MODIFIED="1393943542955" TEXT="Returns the highest mutual information attainable using a grid of x columns and y rows o D"/>
<node CREATED="1393943552046" ID="ID_1547290520" MODIFIED="1393943777683" TEXT="ApproxMaxMI">
<font BOLD="true" NAME="SansSerif" SIZE="12"/>
<node CREATED="1393943565724" ID="ID_1083616701" MODIFIED="1393943773530" TEXT="OptimizeXAxis">
<font BOLD="true" NAME="SansSerif" SIZE="12"/>
<node CREATED="1393943584057" ID="ID_902526452" MODIFIED="1393943610331" TEXT="Given a fixed y-axis partition, OptimizeXAxis finds the x-axis partition that will lead to a grid with maximal mutual information"/>
<node CREATED="1393943677006" HGAP="41" ID="ID_202457643" MODIFIED="1393944111509" TEXT="Draws x-axis partition lines only between runs of consecutive points that all in the same row of the y-axis partition (clumps)" VSHIFT="38"/>
<node CREATED="1393943710069" HGAP="25" ID="ID_1030166625" MODIFIED="1393944107749" TEXT="GetClumpsPartition" VSHIFT="47">
<font BOLD="true" NAME="SansSerif" SIZE="12"/>
<node CREATED="1393943718105" ID="ID_153948474" MODIFIED="1393943766527" TEXT="finds the Edged of clumps: it returns he minimal partition that separates every pair of points that lie ind distinct clumps"/>
</node>
<node CREATED="1393944000818" HGAP="27" ID="ID_372319453" MODIFIED="1393944170694" TEXT="GetSuperClumpsPartition" VSHIFT="73">
<font BOLD="true" NAME="SansSerif" SIZE="12"/>
<node CREATED="1393944026136" ID="ID_1804122753" MODIFIED="1393944057823" TEXT="used to reduce the runtime complexity of OptimizeXAxis"/>
<node CREATED="1393944058707" ID="ID_274139596" MODIFIED="1393944250173" TEXT="merges clumps in superclumps in a way that aims to have each superclump contain approximately the same number of points. The algorithm then forgets about the clumps and only considers drawing grid-lines between the supeclumps."/>
<node CREATED="1393944170695" ID="ID_484624288" MODIFIED="1393944284280" TEXT="the parameter k_hat allows for a standard efficiency vs. optimality tradeoff."/>
</node>
</node>
</node>
</node>
</node>
</node>
</map>
