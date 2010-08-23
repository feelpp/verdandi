<?php echo '<?xml version="1.0"  encoding="iso-8859-1"?'.'>' ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<?php $root='..';?>

<head>
<meta http-equiv="Content-Type" content="text/html;charset=UTF-8">
<title>Verdandi user's guide</title>
<link rel="stylesheet" type="text/css" href="<?php echo $root?>/content.css">
<link rel="stylesheet" href="tabs.css" type="text/css">
<link rel="stylesheet" href="guide.css" type="text/css">
<script type="text/javascript" src="prettify.js"></script>
</head>

<body onload="prettyPrint()">

<div class="page">

<?php #include $root.'/header.php';?>

<div class="doc">

<?php function HL($file_, $section_, $string_)
{
if ($file_ == $section_)
  echo '<em>'.$string_.' </em>';
else
  echo '<a href="'.$section_.'.php">'.$string_.'</a>';
}; ?>

<?php $file=basename($_SERVER['REQUEST_URI'], ".php"); $file = explode(".", $file); $file = $file[0];?>

<div class="nav">

<ul>
<li class="jelly"> <b>USER'S GUIDE</b> </li>
<li class="jelly"> <?php HL($file, "index", "Introduction");?>  </li>

<li class="jelly"> <?php HL($file, "getting_started", "Getting Started");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "getting_started"
or basename($_SERVER['REQUEST_URI'], ".php") == "installation"
or basename($_SERVER['REQUEST_URI'], ".php") == "notation"
or basename($_SERVER['REQUEST_URI'], ".php") == "overview")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "installation", "Installation");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "notation", "Notation");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "overview", "Overview");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "assimilation_methods", "Assimilation Methods");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "assimilation_methods"
or basename($_SERVER['REQUEST_URI'], ".php") == "optimal_interpolation")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "optimal_interpolation", "Optimal Interpolation");
  echo '</li> </ul>';
} ?>

</li>


<li class="jelly"> <?php HL($file, "models", "Models");?> </li>

<li class="jelly"> <?php HL($file, "observations", "Observations");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "observations"
or basename($_SERVER['REQUEST_URI'], ".php") == "linear_observation_manager")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "linear_observation_manager", "Linear Observation Manager");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "plugging_in_verdandi", "Plugging in Verdandi");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "plugging_in_verdandi"
or basename($_SERVER['REQUEST_URI'], ".php") == "plugging_model"
or basename($_SERVER['REQUEST_URI'], ".php") == "plugging_observation")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "plugging_model", "Model");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "plugging_observation", "Observations");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "tips", "Tips");?>  </li>

<li class="jelly"> <?php HL($file, "debugging", "Debugging");?>  </li>

<li class="jelly"> <?php HL($file, "tutorial", "Tutorial");?>  </li>

<li class="jelly"> <?php HL($file, "python", "Python");?> </li>
<li class="jelly"> <b>API REFERENCE</b> </li>
<li class="jelly"> <?php HL($file, "annotated", "Classes");?>
<ul class="navsubul"> <li class="jelly"> <?php HL($file, "annotated", "Class List");?> </li> 
<li class="jelly"> <?php HL($file, "hierarchy", "Class Hierarchy");?> </li>
<li class="jelly"> <?php HL($file, "functions", "Class Members");?>
</li> </ul> </li>
<li class="jelly"> <?php HL($file, "namespacemembers", "Functions");?> </li>
<li class="jelly"> Search for <form action="search.php" method="get">
    <input class="search" type="text" name="query" value="" size="20" accesskey="s">
  </form>
</li>
<!-- <li class="jelly"> <?php HL($file, "faq", "F.A.Q.");?> </li>-->
<li class="jelly"> <a
href="mailto:seldon-help@lists.sourceforge.net"
style="color:black">Support</a></li>
</ul>

</div>

<div class="doxygen">
