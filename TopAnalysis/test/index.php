<html>
<head>
    <title><?php echo getcwd(); ?></title>
    <style type='text/css'>
        body {
            font-family: "Candara", sans-serif;
            font-size: 9pt;
            line-height: 10.5pt;
        }
        div.pic h3 {
            font-size: 11pt;
            margin: 0.5em 1em 0.2em 1em;
        }
        div.pic p {
            font-size: 11pt;
            margin: 0.2em 1em 0.1em 1em;
        }
        div.pic {
            display: block;
            float: left;
            background-color: white;
            border: 1px solid #ccc;
            padding: 2px;
            text-align: center;
            margin: 2px 10px 10px 2px;
            -moz-box-shadow: 7px 5px 5px rgb(80,80,80);    /* Firefox 3.5 */
            -webkit-box-shadow: 7px 5px 5px rgb(80,80,80); /* Chrome, Safari */
            box-shadow: 7px 5px 5px rgb(80,80,80);         /* New browsers */
        }
        a { text-decoration: none; color: rgb(80,0,0); }
        a:hover { text-decoration: underline; color: rgb(255,80,80); }
        div.dirlinks h2 {  margin-bottom: 4pt; margin-left: -24pt; color: rgb(80,0,0);  }
        div.dirlinks {  margin: 0 24pt; }
        div.dirlinks a {
            font-size: 11pt; font-weight: bold;
            padding: 0 0.5em;
        }

        .categorylink{
            font-size: 20pt;
            font-weight: bold;
            height: 50px;
            display: inline-block;
            margin-top: 30px;
        }
        .categorydiv{
            float: left;
        }

    </style>
<script language="JavaScript">
    function toggle(id) {
        var state = document.getElementById(id).style.display;
            if (state == 'block') {
                document.getElementById(id).style.display = 'none';
            } else {
                document.getElementById(id).style.display = 'block';
        }
    }
</script>
</head>
<body>
    <h1><?php echo getcwd(); ?></h1>
    <?php
    $has_subs = true;
    foreach (glob("*") as $filename) {
        if (is_dir($filename) && !preg_match("/^\..*|.*private.*/", $filename)) {
            $has_subs = true;
            break;
        }
    }
    if ($has_subs) {
        print "<div class=\"dirlinks\">\n";
        print "<h2>Directories</h2>\n";
        print "<a href=\"../\">[parent]</a> ";
        foreach (glob("*") as $filename) {
            if (is_dir($filename) && !preg_match("/^\..*|.*private.*/", $filename)) {
                print " <a href=\"$filename\">[$filename]</a>";
            }
        }
        print "</div>";
    }

    foreach (array("00_README.txt", "README.txt", "readme.txt") as $readme) {
        if (file_exists($readme)) {
            print "<pre class='readme'>\n"; readfile($readme); print "</pre>";
        }
    }
    ?>

    <h2><a name="plots">Plots</a></h2>
    <p><form>Filter: <input type="text" name="match" size="30" value="<?php if (isset($_GET['match'])) print htmlspecialchars($_GET['match']);  ?>" /><input type="Submit" value="Go" /></form></p>
    <div>
         <?php
         $categories = array();
         $othercategories = array();
         $filenames = glob("*.png"); sort($filenames);
         foreach ($filenames as $filename) {
             $parts = explode("_", $filename);
             $category = $parts[0];
             if( !in_array( $category , $categories ) ){
                 if( isset( $_GET['cats'] ) ){
                     $cats = explode( ":" , $_GET['cats'] );
                     if( in_array( $category , $cats ) ){
                         array_push( $categories , $category );
                     }else if(!in_array( $category , $othercategories )){
                         array_push( $othercategories , $category );
                     }
                 }else if( isset( $_GET['catfilter'] ) ){
                     if( fnmatch('*'.$_GET['catfilter'].'*', $category ) ){
                         array_push( $categories , $category );
                     }else if(!in_array( $category , $othercategories )){
                         array_push( $othercategories , $category );
                     }
                 }
                 else{
                     array_push( $categories , $category );
                 }
             }
         }
        $displayed = array();
        if ($_GET['noplots']) {
            print "Plots will not be displayed.\n";
        } else {
            $other_exts = array('.pdf', '.cxx', '.eps', '.root', '.txt');
            foreach ($categories as $category) {
                print "<div><a href=\"#\" class=\"categorylink\" onclick=\"toggle('div$category')\">Category: $category</a></br>";
                print "<div id='div$category' class='categorydiv' style=\"display:none;\">";
            $filenames = glob("*.png"); sort($filenames);
            foreach ($filenames as $filename) {
                $parts = explode("_", $filename);
                $category_ = $parts[0];
                if ( $category != $category_ ) continue;
                if (isset($_GET['match']) && !fnmatch('*'.$_GET['match'].'*', $filename)) continue;
                array_push($displayed, $filename);
                print "<div class='pic'>\n";
                print "<h3><a href=\"$filename\">$filename</a></h3>";
                print "<a href=\"$filename\"><img src=\"$filename\" style=\"border: none; width: 400px; \"></a>";
                $others = array();
                foreach ($other_exts as $ex) {
                    $other_filename = str_replace('.png', $ex, $filename);
                    if (file_exists($other_filename)) {
                        array_push($others, "<a class=\"file\" href=\"$other_filename\">[" . $ex . "]</a>");
                        if ($ex != '.txt') array_push($displayed, $other_filename);
                    }
                }
                if ($others) print "<p>Also as ".implode(', ',$others)."</p>";
                print "</div>";
            }
            print "</div></div>";
            }
        }
          ?>
    </div>
    <div style="display: block; clear:both;">
        <h2><a name="files">Other files/categories</a></h2>
        <ul>
            <?
             if( isset( $_GET['cats'] ) || isset( $_GET['catfilter'] ) ){
                 foreach( $othercategories as $cat ){
                     print "<li>[CAT] <a href=?cats=$cat>$cat</a></li>";
                 }
             }else{
            foreach (glob("*") as $filename) {
                if ($_GET['noplots'] || !in_array($filename, $displayed)) {
                    if (isset($_GET['match']) && !fnmatch('*'.$_GET['match'].'*', $filename)) continue;
                    if ($filename === 'index.php') continue;
                    if (is_dir($filename)) {
                        print "<li>[DIR] <a href=\"$filename\">$filename</a></li>";
                    } else {
                        print "<li><a href=\"$filename\">$filename</a></li>";
                    }
                }
            }
             }
            ?>
        </ul>
    </div>
</body>
</html>
