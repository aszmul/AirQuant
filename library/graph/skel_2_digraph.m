function [digraphout, glinkout, gnode, isloops] = skel_2_digraph(skel, method)
    % Generate digraph from skeleton with loop breaking
    if nargin < 2
        method = 'topnode';
    end
    
    isloops = 0;

    % Get initial undirected graph
    [gadj, gnode, glink] = Skel2Graph3D(skel, 1);
    
    % Create initial graph and check connectivity
    G_initial = digraph(gadj);
    [comp_ids_initial, comp_sizes_initial] = conncomp(G_initial, 'Type', 'weak');
    fprintf('Initial state: %d connected components\n', length(unique(comp_ids_initial)));
    
    % Check for loops
    G_loop_check = digraph(gadj);
    if length(glink) ~= height(G_loop_check.Edges)/2
        warning('Loops detected in skeleton. Will break loops while preserving connectivity.')
        isloops = 1;
        
        % Find and break loops using modified minimum spanning tree
        G_undirected = graph(gadj);
        
        % Use edge lengths as weights (longer edges more likely to be preserved)
        weights = ones(size(G_undirected.Edges.EndNodes, 1), 1);
        for i = 1:length(glink)
            weight_idx = findedge(G_undirected, glink(i).n1, glink(i).n2);
            if weight_idx > 0
                weights(weight_idx) = 1/length(glink(i).point); % Inverse length to prefer keeping shorter edges
            end
        end
        G_undirected.Edges.Weight = weights;
        
        % Get initial component info
        [bins_before, ~] = conncomp(G_undirected);
        num_comp_before = length(unique(bins_before));
        
        % Get minimum spanning tree while tracking removed edges
        T = minspantree(G_undirected, 'Method', 'sparse');
        edges_to_remove = [];
        
        % Carefully remove edges that create loops
        for i = 1:length(glink)
            if ~findedge(T, glink(i).n1, glink(i).n2)
                % Test if removing this edge would increase number of components
                temp_G = G_undirected;
                temp_G = rmedge(temp_G, findedge(temp_G, glink(i).n1, glink(i).n2));
                [bins_after, ~] = conncomp(temp_G);
                num_comp_after = length(unique(bins_after));
                
                if num_comp_after <= num_comp_before
                    edges_to_remove = [edges_to_remove; i];
                end
            end
        end
        
        % Remove identified edges
        fprintf('Removing %d edges to break loops\n', length(edges_to_remove));
        glink(edges_to_remove) = [];
        
        % Check connectivity after loop breaking
        edges_remaining = [[glink.n1]', [glink.n2]'];
        G_after = graph(edges_remaining(:,1), edges_remaining(:,2));
        [comp_ids_after, comp_sizes_after] = conncomp(G_after);
        fprintf('After loop breaking: %d connected components\n', length(unique(comp_ids_after)));
        
        if length(unique(comp_ids_after)) > length(unique(comp_ids_initial))
            warning('Loop breaking increased number of components! This needs manual review.');
            fprintf('Component size distribution before:\n');
            disp(tabulate(comp_ids_initial));
            fprintf('Component size distribution after:\n');
            disp(tabulate(comp_ids_after));
        end
    end

    % Create bidirectional edges
    edges_og = [[glink.n1]', [glink.n2]'];
    edges_rev = [[glink.n2]', [glink.n1]'];
    weights = cellfun(@numel, {glink.point});
    weights = [weights'; weights'];
    edges_twoway = [edges_og; edges_rev];
    Edgetable = table(edges_twoway, weights, 'VariableNames', {'EndNodes', 'Weight'});
    G = digraph(Edgetable);
    nodelist = 1:numnodes(G);

    % Find connected components and origin nodes
    [bins, binsize] = conncomp(G, 'Type', 'weak');
    originnode = zeros(max(bins), 1);

    if strcmp(method, 'topnode') 
        for ii = 1:max(bins)
            binbool = (bins==ii);
            binidx = nodelist(binbool);
            gnodeii = gnode(binidx);
            [~, binorigin] = max([gnodeii.comz]);
            originnode(ii) = binidx(binorigin);
        end
    elseif isnumeric(method)
        assert(length(unique(bins)) < 2, 'Can only be used with one connected component.')
        assert(length(method)==3, 'Must be of length 3')
        node_coords = [[gnode.comy]; [gnode.comx]; [gnode.comz]]';
        originnode = dsearchn(node_coords, method);
    else
        error('Chosen method is invalid.')
    end

    % BF search from origin node
    allnode_discovery = cell(max(bins), 1);
    for ii = 1:max(bins)
        allnode_discovery{ii} = bfsearch(G, originnode(ii));
    end
    node_discovery = cell2mat(allnode_discovery);
    
    % Reorder nodes
    G = reordernodes(G, node_discovery);
    gnode = gnode(node_discovery);
    
    % Update node connections
    for iinode = 1:length(gnode)
        conns = gnode(iinode).conn;
        for iiconn = 1:length(conns)
            gnode(iinode).conn(iiconn) = find(node_discovery == conns(iiconn));
        end
    end
    
    % Update link endpoints
    for iilink = 1:length(glink)
        glink(iilink).n1 = find(node_discovery == glink(iilink).n1);
        glink(iilink).n2 = find(node_discovery == glink(iilink).n2);
    end

    % Create outward facing digraph
    removal = zeros(height(G.Edges)/2, 1);
    j = 1;
    for i = 1:height(G.Edges)
        if (G.Edges.EndNodes(i,1) - G.Edges.EndNodes(i,2)) > 0
            removal(j) = i;
            j = j + 1;
        end
    end

    % Remove edges in opposite direction
    removal(removal == 0) = [];
    G = rmedge(G, removal);

    % Ensure glink directions match digraph
    for i = 1:length(glink)
        glink_nodes = [glink(i).n1, glink(i).n2];
        if ~ismember(glink_nodes, G.Edges.EndNodes, 'rows')
            glink(i).n1 = glink_nodes(2);
            glink(i).n2 = glink_nodes(1);
            glink(i).point = fliplr(glink(i).point);
        end
    end    

    % Create final digraph
    edges = [[glink.n1]', [glink.n2]'];
    weights = zeros(length(glink), 1);
    for i = 1:length(glink)
        weights(i) = length(glink(i).point);
    end
    labels = (1:length(glink))';
    Edgetable = table(edges, weights, labels, 'VariableNames', {'EndNodes', 'Weight', 'Label'});
    
    digraphout = digraph(Edgetable);
    
    % Create node table with properties
    n_nodes = numnodes(digraphout);
    fprintf('Creating node properties table for %d nodes\n', n_nodes);
    
    % Verify node count matches
    if length(gnode) ~= n_nodes
        warning('Node count mismatch! Graph has %d nodes but gnode has %d entries', ...
            n_nodes, length(gnode));
        % Ensure gnode matches graph nodes
        gnode = gnode(1:n_nodes);
    end
    
    % Create node property arrays
    comx = zeros(n_nodes, 1);
    comy = zeros(n_nodes, 1);
    comz = zeros(n_nodes, 1);
    ep = zeros(n_nodes, 1);
    label = (1:n_nodes)';
    
    % Fill arrays from gnode
    for i = 1:n_nodes
        comx(i) = gnode(i).comx;
        comy(i) = gnode(i).comy;
        comz(i) = gnode(i).comz;
        ep(i) = gnode(i).ep;
    end
    
    % Assign to graph
    digraphout.Nodes.comx = comx;
    digraphout.Nodes.comy = comy;
    digraphout.Nodes.comz = comz;
    digraphout.Nodes.ep = ep;
    digraphout.Nodes.label = label;

    % Process final edge ordering
    try
        % Get components and sort by size
        [comp_ids_final, comp_sizes] = conncomp(digraphout, 'Type', 'weak');
        unique_comps = unique(comp_ids_final);
        [~, sort_idx] = sort(comp_sizes, 'descend');
        unique_comps = unique_comps(sort_idx);
        
        fprintf('Final digraph before cleanup: %d connected components\n', length(unique_comps));
        
        % Initialize tracking of reached edges
        n_edges = size(digraphout.Edges, 1);
        reached_edges = false(n_edges, 1);
        ordered_edges = [];
        
        % Process each component
        for comp = unique_comps'
            comp_nodes = find(comp_ids_final == comp);
            
            % Find root (highest z-coord)
            [~, max_z_idx] = max(digraphout.Nodes.comz(comp_nodes));
            root = comp_nodes(max_z_idx);
            
            % Get edges for this component using BFS
            [~, comp_edges] = bfsearch(digraphout, root, 'edgetonew');
            
            % Mark these edges as reached
            reached_edges(comp_edges) = true;
            ordered_edges = [ordered_edges; comp_edges];
        end
        
        % Check for unreached edges
        unreached = find(~reached_edges);
        if ~isempty(unreached)
            fprintf('Removing %d unreached edges\n', length(unreached));
            
            % Only keep edges that were reached
            if max(ordered_edges) <= length(digraphout.Edges.Label)
                edge_labels = digraphout.Edges.Label(ordered_edges);
                if max(edge_labels) <= length(glink)
                    glinkout = glink(edge_labels);
                else
                    warning('Edge labels exceed glink length. Trimming edges.');
                    valid_idx = edge_labels <= length(glink);
                    glinkout = glink(edge_labels(valid_idx));
                end
            else
                warning('Edge indices exceed graph size. Using valid edges only.');
                valid_idx = ordered_edges <= length(digraphout.Edges.Label);
                edge_labels = digraphout.Edges.Label(ordered_edges(valid_idx));
                glinkout = glink(edge_labels);
            end
            
            % Create new digraph with only reached edges
            edges = [[glinkout.n1]', [glinkout.n2]'];
            weights = zeros(length(glinkout), 1);
            for i = 1:length(glinkout)
                weights(i) = length(glinkout(i).point);
            end
            labels = (1:length(glinkout))';
            Edgetable = table(edges, weights, labels, 'VariableNames', {'EndNodes', 'Weight', 'Label'});
            
            % Create node properties table
            n_nodes = max(max(edges));
            NodeTable = table((1:n_nodes)', ...
                            zeros(n_nodes,1), ...
                            zeros(n_nodes,1), ...
                            zeros(n_nodes,1), ...
                            zeros(n_nodes,1), ...
                            (1:n_nodes)', ...
                            'VariableNames', {'Node', 'comx', 'comy', 'comz', 'ep', 'label'});
            
            % Fill node properties
            for i = 1:n_nodes
                if i <= length(gnode)
                    NodeTable.comx(i) = gnode(i).comx;
                    NodeTable.comy(i) = gnode(i).comy;
                    NodeTable.comz(i) = gnode(i).comz;
                    NodeTable.ep(i) = gnode(i).ep;
                end
            end
            
            % Create new digraph with node and edge properties
            digraphout = digraph(Edgetable, NodeTable);
            
            % Verify final graph properties
            [comp_ids_final, comp_sizes] = conncomp(digraphout, 'Type', 'weak');
            unique_comps = unique(comp_ids_final);
            
            fprintf('\nFinal validation:\n');
            fprintf('- Connected components: %d\n', length(unique_comps));
            fprintf('- Nodes: %d\n', numnodes(digraphout));
            fprintf('- Edges: %d\n', size(digraphout.Edges, 1));
            
            % Check edge consistency
            n1_valid = all(digraphout.Edges.EndNodes(:,1) <= numnodes(digraphout));
            n2_valid = all(digraphout.Edges.EndNodes(:,2) <= numnodes(digraphout));
            if ~n1_valid || ~n2_valid
                warning('Invalid node references in edges!');
            else
                fprintf('- All edge endpoints are valid\n');
            end
            
            % Check if all nodes are reachable from their component root
            unreachable_nodes = false;
            for comp = unique_comps'
                comp_nodes = find(comp_ids_final == comp);
                [~, max_z_idx] = max(digraphout.Nodes.comz(comp_nodes));
                root = comp_nodes(max_z_idx);
                reached = dfsearch(digraphout, root);
                if length(reached) ~= length(comp_nodes)
                    unreachable_nodes = true;
                    break;
                end
            end
            
            if unreachable_nodes
                warning('Some nodes are unreachable from their component root!');
            else
                fprintf('- All nodes are reachable from component roots\n');
            end
            
        else
            glinkout = glink(digraphout.Edges.Label(ordered_edges));
            fprintf('\nNo unreached edges to remove. Graph is clean.\n');
        end
        
    catch ME
        warning('Edge processing failed: %s', ME.message);
        % In case of error, attempt to salvage valid edges
        try
            valid_idx = ordered_edges <= length(digraphout.Edges.Label);
            edge_labels = digraphout.Edges.Label(ordered_edges(valid_idx));
            valid_labels = edge_labels <= length(glink);
            glinkout = glink(edge_labels(valid_labels));
        catch
            warning('Failed to process edges. Graph may be incomplete.');
            glinkout = glink;
        end
    end
end
