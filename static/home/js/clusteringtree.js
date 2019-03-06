function renderTree(x) {
    var r = 1200 / 2;
    var innerRadius = r - 170 // change inner radius of tree with this argument
    var names = 0; // indexing for all nodes


    var cluster = d3.layout.cluster()
        .size([360, innerRadius])
        .sort(function(a, b) {
            return (a.value - b.value) || d3.ascending(a.length, b.length);
        })
        .value(function(d) {
            return d.length;
        })
        .children(function(d) {
            return d.branchset;
        })
        .separation(function(a, b) {
            return 1;
        });

    function project(d) {
        var r = d.y,
            a = (d.x - 90) / 180 * Math.PI;
        return [r * Math.cos(a), r * Math.sin(a)];
    }

    function cross(a, b) {
        return a[0] * b[1] - a[1] * b[0];
    }

    function dot(a, b) {
        return a[0] * b[0] + a[1] * b[1];
    }

    function step(d) {
        var s = project(d.source),
            m = project({
                x: d.target.x,
                y: d.source.y
            }),
            t = project(d.target),
            r = d.source.y,
            sweep = d.target.x > d.source.x ? 1 : 0;
        return (
            "M" + s[0] + "," + s[1] +
            "A" + r + "," + r + " 0 0," + sweep + " " + m[0] + "," + m[1] +
            "L" + t[0] + "," + t[1]);
    }

    var rect_dim = (r - 90) * 2;
    var svg_dim = (r - 90) * 2;

    var wrap = d3.select("#clustering-tree").append("svg")
        .attr("width", svg_dim)
        .attr("height", svg_dim)
        .style("-webkit-backface-visibility", "hidden");

    // Catch mouse events in Safari.
    wrap.append("rect")
        .attr("width", rect_dim)
        .attr("height", rect_dim)
        .attr("fill", "none");

    var vis = wrap.append("g")
        .attr("transform", "translate(" + rect_dim / 2 + "," + rect_dim / 2 + ")");

    var start = null,
        rotate = 0,
        div = document.getElementById("clustering-tree");

    //function to catch XY coordinates for mouse
    /*function mouse(e) {

        return [
            e.pageX - div.offsetLeft - r,
            e.pageY - div.offsetTop - r
        ];
    }

    wrap.on("mousedown", function() {
        wrap.style("cursor", "move");
        start = mouse(d3.event);
        d3.event.preventDefault();
    });
    d3.select(window)
        .on("mouseup", function() {
            if (start) {
                wrap.style("cursor", "auto");
                var m = mouse(d3.event);
                var delta = Math.atan2(cross(start, m), dot(start, m)) * 180 / Math.PI;
                rotate += delta;
                if (rotate > 360) rotate %= 360;
                else if (rotate < 0) rotate = (360 + rotate) % 360;
                start = null;
                wrap.style("-webkit-transform", null);
                vis
                    .attr("transform", "translate(" + rect_dim / 2 + "," + rect_dim / 2 + ")rotate(" + rotate + ")")
                    .selectAll("text")
                    .attr("text-anchor", function(d) {
                        return (d.x + rotate) % 360 < 180 ? "start" : "end";
                    })
                    .attr("transform", function(d) {
                        return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 2) + ")rotate(" + ((d.x + rotate) % 360 < 180 ? 0 : 180) + ")";
                    });
            }
        })
        .on("mousemove", function() {
            if (start) {
                var m = mouse(d3.event);
                var delta = Math.atan2(cross(start, m), dot(start, m)) * 180 / Math.PI;
                wrap.style("-webkit-transform", "rotateZ(" + delta + "deg)");
            }
        });*/

    //Branch length function
    function phylo(n, offset) {
        if (n.length != null) offset += n.length * 115;
        n.y = offset;
        if (n.children)
            n.children.forEach(function(n) {
                phylo(n, offset);
            });
    }

    x = newick.parse(x);
    var nodes = cluster.nodes(x);

    nodes.forEach(function(n) {

        if (n.name == "") {
            n.name = names.toString();
            names++;
        }

    });

    //Uncomment the line below to show branch length
    // phylo(nodes[0], 0);

    var link = vis.selectAll("path.link")
        .data(cluster.links(nodes))
        .enter().append("path")
        .attr("class", "link")
        .attr("d", step);


    //D3 selection to distinguish between inner and leaf nodes
    var node = vis.selectAll("g.node")
        .data(nodes.filter(function(n) {
            return n.x !== undefined;
        }))
        .enter().append("g")
        .attr("class", function(n) {
            if (n.children) {
                return "inner node";
            } else {
                //return "leaf node";
                return ('X' + n.name);
            }
        })
        .attr("transform", function(d) {
            return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")";
        })

    // Add terminal node with coloring
    node.append("circle")
        .attr("r", 5);

    var innernodes = vis.selectAll('g.inner.node')
        .append("circle")
        .attr("r", 5);

    // Adding annotations
    /*for (var x in selectivityinfo){
        var spacer = 8
        if(selectivityinfo[x].indexOf("Gs family") >= 0){
          var leafwithname = vis.selectAll('g.X'+x)
              .append("circle")
              .attr("r", 3.25)
              .style("fill", "blue")
              .attr("transform", "translate(" + (23 + spacer) + ",0)");
        }

        if(selectivityinfo[x].indexOf("Gi/Go family") >= 0){
          var leafwithname = vis.selectAll('g.X'+x)
              .append("circle")
              .attr("r", 3.25)
              .style("fill", "red")
              .attr("transform", "translate(" + (23  + 2*spacer) + ",0)");
        }

        if(selectivityinfo[x].indexOf("Gq/G11 family") >= 0){
          var leafwithname = vis.selectAll('g.X'+x)
              .append("circle")
              .attr("r", 3.25)
              .style("fill", "black")
              .attr("transform", "translate(" + (23 + 3*spacer) + ",0)");
        }

        if( selectivityinfo[x].indexOf("G12/G13 family") >= 0){
          var leafwithname = vis.selectAll('g.X'+x)
              .append("circle")
              .attr("r", 3.25)
              .style("fill", "green")
              .attr("transform", "translate(" + (23 + 4*spacer) + ",0)");
        }
    }*/


    var label = vis.selectAll("text")
        .data(nodes.filter(function(d) {
            return d.x !== undefined && !d.children;
        }))
        .enter().append("text")
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) {
            return d.x < 180 ? "start" : "end";
        })
        .attr("transform", function(d) {
            return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 10) + ")rotate(" + (d.x < 180 ? 0 : 180) + ")";
        })
        .text(function(d) {
            return d.name.replace(/_/g, ' ');
        });


    //click count to change clade color selection and unhighlight subtree nodes
    var click_count = 0;

    //On click event to get subtree of clicked node
    var firstTreeClick = true;
    node.on("click", function(d) {
        if (firstTreeClick){
          firstTreeClick = false;
          $("#CN-button").removeClass("disabled");
          $("#DN-button").removeClass("disabled");
        }
        doubleclade(d);
    });

    //Highlight clade when mouse is over node
    node.on("mouseover", function(h) {
        traverse(h);
    });

    //Unhighlight clade when mouse is out of node
    node.on("mouseout", function(u) {
        vis.selectAll("path.selected")
            .classed("selected", false);
        vis.selectAll("text").style("font-weight", "normal");
    });

    //iterate over links for (source,target) info to highlight subtree path
    function highlightlink(src, tgt) {
        var link = d3.selectAll("path.link")[0].filter(function(d) {
            var flag = (d3.select(d).data()[0].source.name == src && d3.select(d).data()[0].target.name == tgt);
            return flag;
        });

        d3.selectAll(link)
            .classed("selected", true);
    }

    //iterate over labels for leaf nodes info to highlight labels
    function highlightlabel(tgt) {
        var highlight_labels = vis.selectAll("text")[0].filter(function(t) {
            var th = (d3.select(t).data()[0].name == tgt);
            return th;
        });

        d3.selectAll(highlight_labels)
            .style("font-weight", "bold");
    }

    //Recursive function to highlight all links of a subtree
    function traverse(node) {
        if (node.children) {
            node.children.forEach(function(d) {
                highlightlabel(d.name);
                highlightlink(node.name, d.name);
                traverse(d);
            });
        }
    }

    //Recursive function to traverse a subtree and deselect all nodes
    function deselect(node) {
        node.isSelected = false;
        if (node.children) {
            node.children.forEach(function(d) {
                //Select all children
                d.isSelected = false;
                deselect(d);
            });
        }
    }


    //Recursive function to traverse a node and its subtree
    function selectSubtree(node, subtree) {
        node.isSelected = true;
        if (node.children) {
            node.children.forEach(function(d) {
                d.isSelected = true;
                highlightlabel(d.name)

                if (d.name.length > 3) {
                    subtree.push(d.name);
                }

                selectSubtree(d, subtree);
            });
        }

        return subtree;
    }

    //Double clade selection
    function doubleclade(node) {
        var groupId = click_count % 2;
        var subtreeArray = selectSubtree(node, []);

        // Reset previous selection - remove class
        vis.selectAll('circle.group' + groupId)
            .classed('group' + groupId, false)
        vis.selectAll('path.group' + groupId)
            .classed('group' + groupId, false)

        var subnode = vis.selectAll('g.inner.node')
            .attr("id", function(d) {
                if (d.isSelected == true) {
                    return 'sub' + groupId;
                }
            });

        // Color new selection - add class
        vis.selectAll('#sub' + groupId).selectAll("circle")
            .attr("class", "group" + groupId);
        vis.selectAll('path.selected')
            .attr("class", "link group" + groupId);

        // Deselect
        vis.selectAll('g#sub' + groupId +'.inner.node')
            .attr("id", "deselect");
        vis.selectAll("text").style("font-weight", "normal");

        makeUL(subtreeArray, click_count % 2);
        deselect(node);
        click_count++;
    }

    //Fill text area with subtree data
    function makeUL(descendents, group) {
        var textarea = document.getElementById('input-targets-'+group);
        textarea.value = descendents.join("\n");
    }
}